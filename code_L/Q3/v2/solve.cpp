// Q3 v3 solver (C++). Owner: 罗懿 (L, 建模手).
// 模型依据：模型/Q3/v3/数学模型.md 式(Q3-v3-1)~(Q3-v3-7)（牛顿—欧拉线性模型 L1）
// 控制方程：M·ÿ + C·ẏ + K·y = f(t)，状态顺序 (z_f, z_o, θ_f, θ_o)
// 数据来源：题目信息/A题/附件3.xlsx（ω=1.7152行）、附件4.xlsx（全部数值从文件读取，无手写数据）

#include <cmath>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

namespace {

const double PI = 3.14159265358979323846;
const double OMEGA_TARGET = 1.7152;   // 题面给定入射波圆频率，用于在附件3中定位对应行
const double C_L = 10000.0;           // 题面给定直线 PTO 阻尼系数 (N·s/m)
const double C_THETA = 1000.0;        // 题面给定旋转 PTO 阻尼系数 (N·m·s/rad)
const double GEOM_BOUND = -0.2019575; // 振子—转轴几何约束下界 (m)：z_o-z_f > -0.2019575
const double ANGLE_LIMIT = 0.20;      // 小角度适用域阈值 (rad)
const double SLIP_LIMIT = 0.20;       // 小滑移适用域阈值 |u|/ℓ_e

struct Params {
    // 附件4：结构/几何/弹簧参数
    double m_f, m_o, R_f, H_c, H_n, R_o, H_o;
    double rho, g, k, l0, k_theta, K_theta;
    // 附件3：水动力系数（ω=1.7152 对应行）
    double omega, m_a, J_a, b_h, b_theta, F_amp, L_amp;
    // 派生锁定常数（由 load_params 一次算得后固定）
    double Kh, z_G, d, I_f, I_o, l_e;
};

// ====================================================================
// XLSX 读取器（手写 ZIP 解析 + XML 提取，仅依赖 zlib）
// ====================================================================
std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("无法打开文件: " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::string inflate_raw(const std::string& in) {
    z_stream zs{};
    if (inflateInit2(&zs, -MAX_WBITS) != Z_OK) throw std::runtime_error("inflateInit 失败");
    zs.next_in = (Bytef*)in.data();
    zs.avail_in = (uInt)in.size();
    std::string out;
    char buf[65536];
    int ret;
    do {
        zs.next_out = (Bytef*)buf;
        zs.avail_out = sizeof(buf);
        ret = inflate(&zs, Z_NO_FLUSH);
        out.append(buf, zs.total_out - out.size());
        if (ret == Z_STREAM_END) break;
        if (ret != Z_OK) { inflateEnd(&zs); throw std::runtime_error("inflate 失败"); }
    } while (zs.avail_out == 0);
    inflateEnd(&zs);
    return out;
}

struct Cell { bool is_str = false; double num = 0.0; std::string str; };

std::vector<std::string> parse_shared_strings(const std::string& xml) {
    std::vector<std::string> ss;
    size_t i = 0;
    while ((i = xml.find("<si>", i)) != std::string::npos) {
        size_t end = xml.find("</si>", i);
        if (end == std::string::npos) break;
        std::string si = xml.substr(i + 4, end - (i + 4));
        std::string text;
        size_t j = 0;
        while ((j = si.find("<t", j)) != std::string::npos) {
            size_t ts = si.find(">", j);
            size_t te = si.find("</t>", ts);
            if (ts == std::string::npos || te == std::string::npos) break;
            text += si.substr(ts + 1, te - (ts + 1));
            j = te + 4;
        }
        ss.push_back(text);
        i = end + 5;
    }
    return ss;
}

std::map<std::string, Cell> parse_sheet(const std::string& xml, const std::vector<std::string>& ss) {
    std::map<std::string, Cell> cells;
    size_t i = 0;
    while ((i = xml.find("<c ", i)) != std::string::npos) {
        size_t cend = xml.find("</c>", i);
        if (cend == std::string::npos) break;
        std::string c = xml.substr(i, cend - i + 4);
        size_t rpos = c.find(" r=\"");
        if (rpos == std::string::npos) { i = cend + 4; continue; }
        size_t rs = rpos + 4;
        size_t re = c.find("\"", rs);
        std::string ref = c.substr(rs, re - rs);
        Cell cell;
        cell.is_str = (c.find(" t=\"s\"") != std::string::npos)
                     || (c.find(" t=\"inlineStr\"") != std::string::npos);
        size_t vpos = c.find("<v>");
        if (vpos != std::string::npos) {
            size_t ve = c.find("</v>", vpos);
            std::string v = c.substr(vpos + 3, ve - (vpos + 3));
            if (cell.is_str) {
                int idx = std::stoi(v);
                cell.str = (idx >= 0 && idx < (int)ss.size()) ? ss[idx] : "";
            } else {
                cell.num = std::stod(v);
            }
        } else {
            size_t ip = c.find("<is>");
            if (ip != std::string::npos) {
                size_t ts = c.find("<t", ip);
                size_t tse = c.find(">", ts);
                size_t te = c.find("</t>", tse);
                cell.str = c.substr(tse + 1, te - (tse + 1));
            }
        }
        cells[ref] = cell;
        i = cend + 4;
    }
    return cells;
}

std::map<std::string, std::string> unzip_entries(const std::string& data) {
    std::map<std::string, std::string> out;
    size_t eocd = std::string::npos;
    for (size_t i = data.size() >= 22 ? data.size() - 22 : 0; ; --i) {
        if (i + 4 <= data.size() &&
            (uint8_t)data[i] == 0x50 && (uint8_t)data[i + 1] == 0x4b &&
            (uint8_t)data[i + 2] == 0x05 && (uint8_t)data[i + 3] == 0x06) { eocd = i; break; }
        if (i == 0) break;
    }
    if (eocd == std::string::npos) throw std::runtime_error("不是有效的 xlsx(zip) 文件");
    uint32_t co = (uint32_t)(uint8_t)data[eocd + 16]
        | ((uint32_t)(uint8_t)data[eocd + 17] << 8)
        | ((uint32_t)(uint8_t)data[eocd + 18] << 16)
        | ((uint32_t)(uint8_t)data[eocd + 19] << 24);
    size_t p = co;
    while (p + 4 <= data.size() &&
           (uint8_t)data[p] == 0x50 && (uint8_t)data[p + 1] == 0x4b &&
           (uint8_t)data[p + 2] == 0x01 && (uint8_t)data[p + 3] == 0x02) {
        uint16_t m = (uint16_t)((uint8_t)data[p + 10] | ((uint8_t)data[p + 11] << 8));
        uint32_t cs = (uint32_t)(uint8_t)data[p + 20]
            | ((uint32_t)(uint8_t)data[p + 21] << 8)
            | ((uint32_t)(uint8_t)data[p + 22] << 16)
            | ((uint32_t)(uint8_t)data[p + 23] << 24);
        uint16_t fnl = (uint16_t)((uint8_t)data[p + 28] | ((uint8_t)data[p + 29] << 8));
        uint16_t xl = (uint16_t)((uint8_t)data[p + 30] | ((uint8_t)data[p + 31] << 8));
        uint16_t cl = (uint16_t)((uint8_t)data[p + 32] | ((uint8_t)data[p + 33] << 8));
        uint32_t lo = (uint32_t)(uint8_t)data[p + 42]
            | ((uint32_t)(uint8_t)data[p + 43] << 8)
            | ((uint32_t)(uint8_t)data[p + 44] << 16)
            | ((uint32_t)(uint8_t)data[p + 45] << 24);
        std::string name = data.substr(p + 46, fnl);
        size_t lp = lo;
        uint16_t lfnl = (uint16_t)((uint8_t)data[lp + 26] | ((uint8_t)data[lp + 27] << 8));
        uint16_t lxl = (uint16_t)((uint8_t)data[lp + 28] | ((uint8_t)data[lp + 29] << 8));
        size_t ds = lp + 30 + lfnl + lxl;
        std::string comp = data.substr(ds, cs);
        out[name] = (m == 0) ? comp : inflate_raw(comp);
        p += 46 + fnl + xl + cl;
    }
    return out;
}

std::map<std::string, Cell> read_xlsx(const std::string& path) {
    std::string data = read_file(path);
    auto entries = unzip_entries(data);
    auto sit = entries.find("xl/sharedStrings.xml");
    std::vector<std::string> ss;
    if (sit != entries.end()) ss = parse_shared_strings(sit->second);
    auto it = entries.find("xl/worksheets/sheet1.xml");
    if (it == entries.end()) throw std::runtime_error("xlsx 缺少 sheet1.xml: " + path);
    return parse_sheet(it->second, ss);
}

std::string col_letters(const std::string& ref) {
    size_t i = 0;
    while (i < ref.size() && std::isalpha((unsigned char)ref[i])) ++i;
    return ref.substr(0, i);
}
int row_number(const std::string& ref) {
    size_t i = 0;
    while (i < ref.size() && std::isalpha((unsigned char)ref[i])) ++i;
    return std::stoi(ref.substr(i));
}
std::string norm_label(const std::string& s) {
    size_t pos = s.rfind(" (");
    if (pos != std::string::npos) return s.substr(0, pos);
    return s;
}

Params load_params(const std::string& attach3, const std::string& attach4) {
    // 附件4：A 列为参数名(含单位)，B 列为数值
    auto c4 = read_xlsx(attach4);
    std::map<std::string, double> m4;
    for (auto& kv : c4) {
        if (col_letters(kv.first) == "A" && kv.second.is_str) {
            int r = row_number(kv.first);
            std::string refB = "B" + std::to_string(r);
            auto it = c4.find(refB);
            if (it != c4.end() && !it->second.is_str) m4[norm_label(kv.second.str)] = it->second.num;
        }
    }
    auto need4 = [&](const std::string& name) -> double {
        auto it = m4.find(name);
        if (it == m4.end()) throw std::runtime_error("附件4 缺少参数: " + name);
        return it->second;
    };
    Params p;
    p.m_f = need4("浮子质量");
    p.m_o = need4("振子质量");
    p.R_f = need4("浮子底半径");
    p.R_o = need4("振子半径");
    p.H_o = need4("振子高度");
    p.H_n = need4("浮子圆锥部分高度");
    p.H_c = need4("浮子圆柱部分高度");
    p.rho = need4("海水的密度");
    p.g = need4("重力加速度");
    p.k = need4("弹簧刚度");
    p.l0 = need4("弹簧原长");
    p.k_theta = need4("扭转弹簧刚度");
    p.K_theta = need4("静水恢复力矩系数");

    // 附件3：B 列为频率 ω，定位 ω=1.7152 的行，读取实际频率及水动力系数
    auto c3 = read_xlsx(attach3);
    bool found = false;
    for (auto& kv : c3) {
        if (col_letters(kv.first) == "B" && !kv.second.is_str) {
            int r = row_number(kv.first);
            if (std::fabs(kv.second.num - OMEGA_TARGET) < 5e-4) {
                auto get3 = [&](const std::string& col) -> double {
                    std::string ref = col + std::to_string(r);
                    auto it = c3.find(ref);
                    if (it == c3.end() || it->second.is_str)
                        throw std::runtime_error("附件3 缺少列 " + col);
                    return it->second.num;
                };
                p.omega = get3("B");
                p.m_a = get3("C");
                p.J_a = get3("D");
                p.b_h = get3("E");
                p.b_theta = get3("F");
                p.F_amp = get3("G");
                p.L_amp = get3("H");
                found = true;
                break;
            }
        }
    }
    if (!found) throw std::runtime_error("附件3 未找到 ω=1.7152 对应的 Q3 行");

    // ---- 派生锁定常数，一次算得后固定（与 v3 模型锁定值一致）----
    p.Kh = p.rho * p.g * PI * p.R_f * p.R_f;                       // 垂荡静水恢复刚度

    double S_t = PI * p.R_f * p.R_f;
    double S_c = 2.0 * PI * p.R_f * p.H_c;
    double S_n = PI * p.R_f * std::sqrt(p.R_f * p.R_f + p.H_n * p.H_n);
    double S_total = S_t + S_c + S_n;
    double sigma = p.m_f / S_total;
    double m_t = sigma * S_t;
    double m_c = sigma * S_c;
    double m_n = sigma * S_n;
    double z_t = p.H_n + p.H_c;
    double z_c = p.H_n + p.H_c * 0.5;
    double z_n = (2.0 / 3.0) * p.H_n;
    p.z_G = (m_t * z_t + m_c * z_c + m_n * z_n) / p.m_f;           // 浮子质心高度
    p.d = p.H_n - p.z_G;                                           // 转轴偏置（有符号）
    p.I_f = m_t * (p.R_f * p.R_f / 4.0 + (z_t - p.z_G) * (z_t - p.z_G))
           + m_c * (p.R_f * p.R_f / 2.0 + p.H_c * p.H_c / 12.0 + (z_c - p.z_G) * (z_c - p.z_G))
           + m_n * (p.R_f * p.R_f / 4.0 + p.H_n * p.H_n / 18.0 + (z_n - p.z_G) * (z_n - p.z_G));
    p.I_o = p.m_o / 12.0 * (3.0 * p.R_o * p.R_o + p.H_o * p.H_o);
    p.l_e = p.l0 + p.H_o * 0.5 - p.m_o * p.g / p.k;

    return p;
}

// ====================================================================
// 4x4 线性方程组求解（带列主元高斯消元），用于求 M 的逆
// ====================================================================
bool solve4x4(const double A[4][4], const double b[4], double x[4]) {
    double M[4][4], rhs[4];
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = A[i][j]; rhs[i] = b[i]; }

    for (int col = 0; col < 4; ++col) {
        int best = col;
        for (int row = col + 1; row < 4; ++row)
            if (std::fabs(M[row][col]) > std::fabs(M[best][col])) best = row;
        if (std::fabs(M[best][col]) < 1e-14) return false;
        if (best != col) {
            for (int j = 0; j < 4; ++j) std::swap(M[col][j], M[best][j]);
            std::swap(rhs[col], rhs[best]);
        }
        double piv = M[col][col];
        for (int row = col + 1; row < 4; ++row) {
            double factor = M[row][col] / piv;
            for (int j = col; j < 4; ++j) M[row][j] -= factor * M[col][j];
            rhs[row] -= factor * rhs[col];
        }
    }
    for (int i = 3; i >= 0; --i) {
        double sum = rhs[i];
        for (int j = i + 1; j < 4; ++j) sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }
    return true;
}

void invert4x4(const double A[4][4], double Ainv[4][4]) {
    for (int col = 0; col < 4; ++col) {
        double e[4] = {0, 0, 0, 0};
        e[col] = 1.0;
        double x[4];
        if (!solve4x4(A, e, x))
            throw std::runtime_error("质量矩阵奇异，无法求逆");
        for (int row = 0; row < 4; ++row) Ainv[row][col] = x[row];
    }
}

// ====================================================================
// 2x2 实对称矩阵 [[a,b],[b,d]] 的最小特征值（用于矩阵性质核验）
// ====================================================================
double eig_min_2x2(double a, double b, double d) {
    double tr = a + d;
    double halfdiff = 0.5 * (a - d);
    double disc = std::sqrt(halfdiff * halfdiff + b * b);
    return 0.5 * tr - disc;
}

// ====================================================================
// 常系数块矩阵组装（状态顺序 z_f, z_o, θ_f, θ_o）
// 式 Q3-v3-5：M·ÿ + C·ẏ + K·y = f(t)
// ====================================================================
void build_matrices(const Params& p,
                    double M[4][4], double C[4][4], double K[4][4]) {
    double Mf = p.m_f + p.m_a;    // 浮子等效垂荡质量
    double Jf = p.I_f + p.J_a;    // 浮子等效纵摇惯量
    double d2 = p.d * p.d;
    double le2 = p.l_e * p.l_e;
    double dle = p.d * p.l_e;

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            M[i][j] = C[i][j] = K[i][j] = 0.0;

    // 垂荡块 M_h (Q3-v3-1)
    M[0][0] = Mf;
    M[1][1] = p.m_o;

    // 纵摇块 M_p (Q3-v3-4)
    M[2][2] = Jf + p.m_o * d2;
    M[2][3] = p.m_o * dle;
    M[3][2] = p.m_o * dle;
    M[3][3] = p.I_o + p.m_o * le2;

    // 垂荡块 C_h
    C[0][0] = p.b_h + C_L;
    C[0][1] = -C_L;
    C[1][0] = -C_L;
    C[1][1] = C_L;

    // 纵摇块 C_p
    C[2][2] = p.b_theta + C_THETA;
    C[2][3] = -C_THETA;
    C[3][2] = -C_THETA;
    C[3][3] = C_THETA;

    // 垂荡块 K_h^(2)
    K[0][0] = p.Kh + p.k;
    K[0][1] = -p.k;
    K[1][0] = -p.k;
    K[1][1] = p.k;

    // 纵摇块 K_p
    K[2][2] = p.K_theta + p.k_theta - p.m_o * p.g * p.d;
    K[2][3] = -p.k_theta;
    K[3][2] = -p.k_theta;
    K[3][3] = p.k_theta - p.m_o * p.g * p.l_e;
}

// ====================================================================
// 一阶状态右端函数（8 分量）：y=(z_f,z_o,θ_f,θ_o, ż_f,ż_o,θ̇_f,θ̇_o)
// ====================================================================
void rhs(double t, const double y[8], double dy[8],
         const Params& p,
         const double Minv[4][4], const double C[4][4], const double K[4][4]) {
    const double* q = y;
    const double* qd = y + 4;

    double f[4] = {
        p.F_amp * std::cos(p.omega * t),
        0.0,
        p.L_amp * std::cos(p.omega * t),
        0.0
    };

    double Fv[4];
    for (int i = 0; i < 4; ++i) {
        double sC = 0.0, sK = 0.0;
        for (int j = 0; j < 4; ++j) { sC += C[i][j] * qd[j]; sK += K[i][j] * q[j]; }
        Fv[i] = f[i] - sC - sK;
    }

    double qdd[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            qdd[i] += Minv[i][j] * Fv[j];

    for (int i = 0; i < 4; ++i) dy[i] = qd[i];
    for (int i = 0; i < 4; ++i) dy[4 + i] = qdd[i];
}

void rk4_step(double& t, double y[8], double dt,
              const Params& p,
              const double Minv[4][4], const double C[4][4], const double K[4][4]) {
    double k1[8], k2[8], k3[8], k4[8], yt[8];
    rhs(t, y, k1, p, Minv, C, K);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
    rhs(t + 0.5 * dt, yt, k2, p, Minv, C, K);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
    rhs(t + 0.5 * dt, yt, k3, p, Minv, C, K);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + dt * k3[i];
    rhs(t + dt, yt, k4, p, Minv, C, K);
    for (int i = 0; i < 8; ++i) y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    t += dt;
}

// ====================================================================
// 采样点记录（输出列序与 result3.xlsx 一致：zf,θf,zo,θo）
// ====================================================================
struct Sample {
    double t;        // 时间 (s)
    double zf, vf;   // 浮子垂荡位移 (m)、速度 (m/s)
    double tf, wf;   // 浮子纵摇角 (rad)、角速度 (rad/s)
    double zo, vo;   // 振子绝对垂荡位移 (m)、速度 (m/s)
    double to, wo;   // 振子绝对纵摇角 (rad)、角速度 (rad/s)
    double u, phi;   // 相对位移 u=z_o-z_f (m)、相对角 φ=θ_o-θ_f (rad)
    double Pl, Pth;  // 两类 PTO 瞬时功率 (W)
};

struct FullRec {
    std::vector<Sample> samples;
};

FullRec simulate(const Params& p, double rk_dt = 0.001) {
    double M[4][4], C[4][4], K[4][4];
    build_matrices(p, M, C, K);
    double Minv[4][4];
    invert4x4(M, Minv);

    double T = 2.0 * PI / p.omega;
    double t_max = 40.0 * T;
    double sample_dt = 0.2;
    long sample_interval = (long)std::lround(sample_dt / rk_dt);

    double y[8] = {0};
    double t = 0.0;

    FullRec fr;
    fr.samples.push_back({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    long step = 0;
    while (t < t_max - rk_dt * 0.5) {
        rk4_step(t, y, rk_dt, p, Minv, C, K);
        ++step;
        if (step % sample_interval == 0) {
            double ts = step / (double)sample_interval * sample_dt;
            Sample s;
            s.t = ts;
            s.zf = y[0]; s.vf = y[4];
            s.zo = y[1]; s.vo = y[5];
            s.tf = y[2]; s.wf = y[6];
            s.to = y[3]; s.wo = y[7];
            s.u = y[1] - y[0];
            s.phi = y[3] - y[2];
            double dud = y[5] - y[4];
            double dphid = y[7] - y[6];
            s.Pl = C_L * dud * dud;
            s.Pth = C_THETA * dphid * dphid;
            fr.samples.push_back(s);
        }
    }
    return fr;
}

// ====================================================================
// XLSX 写入器（STORE 方法，无压缩）
// ====================================================================
void put_u16(std::string& s, uint16_t v) {
    s.push_back((char)(v & 0xFF));
    s.push_back((char)((v >> 8) & 0xFF));
}
void put_u32(std::string& s, uint32_t v) {
    s.push_back((char)(v & 0xFF));
    s.push_back((char)((v >> 8) & 0xFF));
    s.push_back((char)((v >> 16) & 0xFF));
    s.push_back((char)((v >> 24) & 0xFF));
}

uint32_t crc32_calc(const std::string& data) {
    uint32_t crc = 0xFFFFFFFF;
    for (unsigned char c : data) {
        crc ^= c;
        for (int k = 0; k < 8; ++k)
            crc = (crc >> 1) ^ (0xEDB88320 & -(crc & 1));
    }
    return ~crc;
}

struct ZipEntry { std::string name, data; };

std::string make_xlsx(const std::vector<ZipEntry>& entries) {
    std::string out, cdir;
    std::vector<uint32_t> offsets;
    for (auto& e : entries) {
        uint32_t crc = crc32_calc(e.data);
        uint32_t usize = (uint32_t)e.data.size();
        offsets.push_back((uint32_t)out.size());

        out += "PK\x03\x04";
        put_u16(out, 20);
        put_u16(out, 0);
        put_u16(out, 0);
        put_u16(out, 0);
        put_u16(out, 0);
        put_u32(out, crc);
        put_u32(out, usize);
        put_u32(out, usize);
        put_u16(out, (uint16_t)e.name.size());
        put_u16(out, 0);
        out += e.name;
        out += e.data;

        cdir += "PK\x01\x02";
        put_u16(cdir, 20);
        put_u16(cdir, 20);
        put_u16(cdir, 0);
        put_u16(cdir, 0);
        put_u16(cdir, 0);
        put_u16(cdir, 0);
        put_u32(cdir, crc);
        put_u32(cdir, usize);
        put_u32(cdir, usize);
        put_u16(cdir, (uint16_t)e.name.size());
        put_u16(cdir, 0);
        put_u16(cdir, 0);
        put_u16(cdir, 0);
        put_u16(cdir, 0);
        put_u32(cdir, 0);
        put_u32(cdir, offsets.back());
        cdir += e.name;
    }

    uint32_t cdir_offset = (uint32_t)out.size();
    out += cdir;

    out += "PK\x05\x06";
    put_u16(out, 0);
    put_u16(out, 0);
    put_u16(out, (uint16_t)entries.size());
    put_u16(out, (uint16_t)entries.size());
    put_u32(out, (uint32_t)cdir.size());
    put_u32(out, cdir_offset);
    put_u16(out, 0);

    return out;
}

void write_result_xlsx(const std::string& path, const std::vector<Sample>& rec) {
    std::string headers[9] = {
        "时间t(s)", "浮子垂荡位移zf(m)", "浮子垂荡速度dzf/dt(m/s)",
        "浮子纵摇角位移θf(rad)", "浮子纵摇角速度dθf/dt(rad/s)",
        "振子垂荡位移zo(m)", "振子垂荡速度dzo/dt(m/s)",
        "振子纵摇角位移θo(rad)", "振子纵摇角速度dθo/dt(rad/s)"
    };
    std::string cols[9] = {"A","B","C","D","E","F","G","H","I"};

    std::ostringstream sheet;
    sheet << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
          << "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
          << "<sheetData>";
    sheet << "<row r=\"1\">";
    for (int j = 0; j < 9; ++j)
        sheet << "<c r=\"" << cols[j] << "1\" t=\"inlineStr\"><is><t>"
              << headers[j] << "</t></is></c>";
    sheet << "</row>";

    for (size_t i = 0; i < rec.size(); ++i) {
        sheet << "<row r=\"" << (i + 2) << "\">";
        sheet << "<c r=\"A" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].t << "</v></c>";
        sheet << "<c r=\"B" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].zf << "</v></c>";
        sheet << "<c r=\"C" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].vf << "</v></c>";
        sheet << "<c r=\"D" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].tf << "</v></c>";
        sheet << "<c r=\"E" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].wf << "</v></c>";
        sheet << "<c r=\"F" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].zo << "</v></c>";
        sheet << "<c r=\"G" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].vo << "</v></c>";
        sheet << "<c r=\"H" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].to << "</v></c>";
        sheet << "<c r=\"I" << (i + 2) << "\"><v>" << std::setprecision(12) << rec[i].wo << "</v></c>";
        sheet << "</row>";
    }
    sheet << "</sheetData></worksheet>";

    std::string ct = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
        "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
        "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
        "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        "</Types>";

    std::string rels = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
        "</Relationships>";

    std::string wb = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">"
        "<sheets><sheet name=\"Sheet1\" sheetId=\"1\" r:id=\"rId1\"/></sheets>"
        "</workbook>";

    std::string wb_rels = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
        "</Relationships>";

    std::vector<ZipEntry> entries;
    entries.push_back({"[Content_Types].xml", ct});
    entries.push_back({"_rels/.rels", rels});
    entries.push_back({"xl/workbook.xml", wb});
    entries.push_back({"xl/_rels/workbook.xml.rels", wb_rels});
    entries.push_back({"xl/worksheets/sheet1.xml", sheet.str()});

    std::string xlsx = make_xlsx(entries);
    std::ofstream f(path, std::ios::binary);
    f.write(xlsx.data(), xlsx.size());
    f.close();
}

void write_result_csv(const std::string& path, const std::vector<Sample>& rec) {
    std::ofstream f(path, std::ios::binary);
    f << "\xEF\xBB\xBF";
    f << "时间t(s),浮子垂荡位移zf(m),浮子垂荡速度dzf/dt(m/s),"
      << "浮子纵摇角位移θf(rad),浮子纵摇角速度dθf/dt(rad/s),"
      << "振子垂荡位移zo(m),振子垂荡速度dzo/dt(m/s),"
      << "振子纵摇角位移θo(rad),振子纵摇角速度dθo/dt(rad/s)\n";
    for (size_t i = 0; i < rec.size(); ++i) {
        f << std::setprecision(12) << rec[i].t << ","
          << rec[i].zf << "," << rec[i].vf << ","
          << rec[i].tf << "," << rec[i].wf << ","
          << rec[i].zo << "," << rec[i].vo << ","
          << rec[i].to << "," << rec[i].wo << "\n";
    }
    f.close();
}

void put_bom(std::ofstream& out) { out << "\xEF\xBB\xBF"; }

// ====================================================================
// 路径解析
// ====================================================================
std::string resolve_data_dir() {
    std::string c[] = {"题目信息/A题", "../../../题目信息/A题", "../题目信息/A题"};
    for (auto& s : c)
        if (std::filesystem::exists(s + "/附件3.xlsx")
            && std::filesystem::exists(s + "/附件4.xlsx")) return s;
    return "题目信息/A题";
}
std::string resolve_out_dir() {
    namespace fs = std::filesystem;
    std::string base = "结果（包括各种过程数据与审查）/运行输出";
    if (!fs::exists(base)) return base + "/run001";
    for (int n = 1; n <= 999; ++n) {
        std::string d = base + "/run" + (n < 10 ? "00" : n < 100 ? "0" : "") + std::to_string(n);
        if (!fs::exists(d)) return d;
    }
    return base + "/run999";
}

int find_idx(double t_target, const std::vector<Sample>& rec) {
    int best = 0;
    double best_diff = 1e30;
    for (int i = 0; i < (int)rec.size(); ++i) {
        double d = std::fabs(rec[i].t - t_target);
        if (d < best_diff) { best = i; best_diff = d; }
    }
    return best;
}

}  // namespace

// ====================================================================
// main
// ====================================================================
int main(int argc, char** argv) {
    std::string data_dir, out_dir;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--data" && i + 1 < argc) data_dir = argv[++i];
        else if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
    }
    if (data_dir.empty()) data_dir = resolve_data_dir();
    if (out_dir.empty()) out_dir = resolve_out_dir();
    std::filesystem::create_directories(out_dir);

    Params p = load_params(data_dir + "/附件3.xlsx", data_dir + "/附件4.xlsx");

    double M[4][4], C[4][4], K[4][4];
    build_matrices(p, M, C, K);
    double Minv[4][4];
    invert4x4(M, Minv);

    std::cout << "======== Q3 v3 参数（牛顿—欧拉线性模型 L1）========" << std::endl;
    std::cout << "ω=" << std::setprecision(10) << p.omega << " s⁻¹  F=" << p.F_amp
              << " N  L=" << p.L_amp << " N·m" << std::endl;
    std::cout << "m_f=" << p.m_f << " m_o=" << p.m_o << " m_a=" << p.m_a
              << " J_a=" << p.J_a << std::endl;
    std::cout << "b_h=" << p.b_h << " b_θ=" << p.b_theta << std::endl;
    std::cout << "K_h=" << p.Kh << " K_θ=" << p.K_theta << std::endl;
    std::cout << "k=" << p.k << " k_θ=" << p.k_theta << std::endl;
    std::cout << "c_l=" << C_L << " c_θ=" << C_THETA << std::endl;
    std::cout << "z_G=" << p.z_G << " d=" << p.d << " ℓ_e=" << p.l_e << std::endl;
    std::cout << "I_f=" << p.I_f << " I_o=" << p.I_o << std::endl;

    double T = 2.0 * PI / p.omega;
    double t_max = 40.0 * T;
    std::cout << "\nT=" << T << " s  t_max=" << t_max << " s (40周期)" << std::endl;

    // ======== 验证指标 1：矩阵对称性与正定/半正定性 ========
    std::cout << "\n====== 矩阵性质核验（验证指标 1）======" << std::endl;
    {
        double max_sym_err = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j) {
                double scale = std::max(std::fabs(M[i][j]), 1.0);
                double e = std::fabs(M[i][j] - M[j][i]) / scale;
                if (e > max_sym_err) max_sym_err = e;
                e = std::fabs(C[i][j] - C[j][i]) / scale;
                if (e > max_sym_err) max_sym_err = e;
                e = std::fabs(K[i][j] - K[j][i]) / scale;
                if (e > max_sym_err) max_sym_err = e;
            }
        std::cout << "最大对称性相对误差 = " << std::setprecision(4) << max_sym_err
                  << " (阈值 1e-10)" << (max_sym_err <= 1e-10 ? "  → 通过" : "  → 不通过") << std::endl;

        double emMh = eig_min_2x2(M[0][0], M[0][1], M[1][1]);
        double emMp = eig_min_2x2(M[2][2], M[2][3], M[3][3]);
        double emCh = eig_min_2x2(C[0][0], C[0][1], C[1][1]);
        double emCp = eig_min_2x2(C[2][2], C[2][3], C[3][3]);
        double emKh = eig_min_2x2(K[0][0], K[0][1], K[1][1]);
        double emKp = eig_min_2x2(K[2][2], K[2][3], K[3][3]);

        std::cout << "min_eig(M_h)=" << std::setprecision(8) << emMh
                  << "  min_eig(M_p)=" << emMp << std::endl;
        std::cout << "min_eig(C_h)=" << emCh << "  min_eig(C_p)=" << emCp << std::endl;
        std::cout << "min_eig(K_h)=" << emKh << "  min_eig(K_p)=" << emKp << std::endl;

        bool pd = (emMh > 0) && (emMp > 0) && (emKh > 0) && (emKp > 0);
        bool psd = (emCh >= 0) && (emCp >= 0);
        std::cout << "正定/半正定判定: " << (pd && psd ? "全部满足" : "存在违反")
                  << (pd ? "" : " [质量/刚度非正定]")
                  << (psd ? "" : " [阻尼非半正定]") << std::endl;

        std::ofstream f(out_dir + "/result3_矩阵特性.csv", std::ios::binary);
        put_bom(f);
        f << "矩阵块,最小特征值,性质要求,判定\n";
        f << "M_h(垂荡质量)," << std::setprecision(12) << emMh << ",正定,"
          << (emMh > 0 ? "PASS" : "FAIL") << "\n";
        f << "M_p(纵摇质量)," << emMp << ",正定,"
          << (emMp > 0 ? "PASS" : "FAIL") << "\n";
        f << "C_h(垂荡阻尼)," << emCh << ",半正定,"
          << (emCh >= 0 ? "PASS" : "FAIL") << "\n";
        f << "C_p(纵摇阻尼)," << emCp << ",半正定,"
          << (emCp >= 0 ? "PASS" : "FAIL") << "\n";
        f << "K_h(垂荡刚度)," << emKh << ",正定,"
          << (emKh > 0 ? "PASS" : "FAIL") << "\n";
        f << "K_p(纵摇刚度)," << emKp << ",正定,"
          << (emKp > 0 ? "PASS" : "FAIL") << "\n";
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_矩阵特性.csv" << std::endl;
    }

    // ======== 验证指标 2：静平衡零状态下初始加速度方向 ========
    {
        double y0[8] = {0};
        double dy0[8];
        rhs(0.0, y0, dy0, p, Minv, C, K);
        std::cout << "\n====== 初始加速度方向核验（验证指标 2）======" << std::endl;
        std::cout << "ÿ(0) = (" << std::setprecision(6)
                  << dy0[4] << ", " << dy0[5] << ", "
                  << dy0[6] << ", " << dy0[7] << ")" << std::endl;
        std::cout << "  浮子垂荡加速度=" << dy0[4] << " m/s² (激励 F>0，向上)"
                  << (dy0[4] > 0 ? " → 方向正确" : " → 方向异常") << std::endl;
        std::cout << "  浮子纵摇角加速度=" << dy0[6] << " rad/s² (激励 L>0)"
                  << (dy0[6] > 0 ? " → 方向正确" : " → 方向异常") << std::endl;
    }

    // ======== 正式积分 ========
    std::cout << "\n======== 正式积分 (L1 线性模型, dt=0.001 s) ========" << std::endl;
    FullRec fr = simulate(p, 0.001);
    std::cout << "输出点数: " << fr.samples.size() << " (含 t=0), 最后时刻 t="
              << fr.samples.back().t << " s" << std::endl;

    // ======== result3.xlsx / result3.csv ========
    write_result_xlsx(out_dir + "/result3.xlsx", fr.samples);
    std::cout << "\n已写入 " << out_dir << "/result3.xlsx" << std::endl;
    write_result_csv(out_dir + "/result3.csv", fr.samples);
    std::cout << "已写入 " << out_dir << "/result3.csv" << std::endl;

    // ======== 论文五时刻 CSV ========
    double key_times[] = {10.0, 20.0, 40.0, 60.0, 100.0};
    {
        std::ofstream f(out_dir + "/result3_key.csv", std::ios::binary);
        put_bom(f);
        f << "时间t(s),浮子垂荡位移zf(m),浮子垂荡速度dzf/dt(m/s),"
          << "浮子纵摇角位移θf(rad),浮子纵摇角速度dθf/dt(rad/s),"
          << "振子垂荡位移zo(m),振子垂荡速度dzo/dt(m/s),"
          << "振子纵摇角位移θo(rad),振子纵摇角速度dθo/dt(rad/s)\n";
        for (double kt : key_times) {
            int idx = find_idx(kt, fr.samples);
            f << std::setprecision(10) << fr.samples[idx].t << ","
              << fr.samples[idx].zf << "," << fr.samples[idx].vf << ","
              << fr.samples[idx].tf << "," << fr.samples[idx].wf << ","
              << fr.samples[idx].zo << "," << fr.samples[idx].vo << ","
              << fr.samples[idx].to << "," << fr.samples[idx].wo << "\n";
        }
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_key.csv" << std::endl;

        std::cout << "\n====== Q3 v3 论文五时刻结果 ======" << std::endl;
        std::cout << "     t      zf     żf     θf    θ̇f    zo      żo      θo    θ̇o" << std::endl;
        for (double kt : key_times) {
            int idx = find_idx(kt, fr.samples);
            std::cout << std::fixed << std::setprecision(6)
                      << std::setw(8) << fr.samples[idx].t
                      << std::setw(9) << fr.samples[idx].zf
                      << std::setw(9) << fr.samples[idx].vf
                      << std::setw(8) << fr.samples[idx].tf
                      << std::setw(8) << fr.samples[idx].wf
                      << std::setw(9) << fr.samples[idx].zo
                      << std::setw(9) << fr.samples[idx].vo
                      << std::setw(8) << fr.samples[idx].to
                      << std::setw(8) << fr.samples[idx].wo << std::endl;
        }
    }

    // ======== 验证指标 6：几何约束 z_o-z_f > -0.2019575 ========
    {
        double min_u = 1e30, max_u = -1e30;
        for (auto& s : fr.samples) {
            if (s.u < min_u) min_u = s.u;
            if (s.u > max_u) max_u = s.u;
        }
        bool ok = (min_u > GEOM_BOUND);
        std::cout << "\n====== 几何约束检查（硬约束 u=z_o-z_f > -0.2019575 m）======" << std::endl;
        std::cout << "min u=" << std::setprecision(10) << min_u << " m  max u=" << max_u << " m" << std::endl;
        std::cout << (ok ? "全程满足几何约束。" : "** 违反几何约束！结果无效！ **") << std::endl;

        std::ofstream f(out_dir + "/result3_几何约束检查.csv", std::ios::binary);
        put_bom(f);
        f << "时间t(s),相对位移u=zo-zf(m),相对角phi=thetao-thetaf(rad),直线PTO功率Pl(W),旋转PTO功率Ptheta(W)\n";
        for (auto& s : fr.samples)
            f << std::setprecision(10) << s.t << "," << s.u << "," << s.phi << ","
              << s.Pl << "," << s.Pth << "\n";
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_几何约束检查.csv (含 u、φ、P_l、P_θ 诊断列)" << std::endl;
    }

    // ======== 验证指标 5：线性化适用域（小角度、小滑移）========
    {
        double max_tf = 0, max_to = 0, max_absu = 0;
        double max_Pl = 0, max_Pth = 0;
        for (auto& s : fr.samples) {
            if (std::fabs(s.tf) > max_tf) max_tf = std::fabs(s.tf);
            if (std::fabs(s.to) > max_to) max_to = std::fabs(s.to);
            if (std::fabs(s.u) > max_absu) max_absu = std::fabs(s.u);
            if (s.Pl > max_Pl) max_Pl = s.Pl;
            if (s.Pth > max_Pth) max_Pth = s.Pth;
        }
        double slip_ratio = max_absu / p.l_e;
        bool angle_ok = (max_tf <= ANGLE_LIMIT) && (max_to <= ANGLE_LIMIT);
        bool slip_ok = (slip_ratio <= SLIP_LIMIT);
        std::cout << "\n====== 线性化适用域核验（验证指标 5）======" << std::endl;
        std::cout << "max|θf|=" << std::setprecision(6) << max_tf << " rad  max|θo|=" << max_to
                  << " rad  (阈值 0.20 rad)  " << (angle_ok ? "满足" : "超限→需升级 L2") << std::endl;
        std::cout << "max|u|/ℓ_e=" << slip_ratio << "  (阈值 0.20)  "
                  << (slip_ok ? "满足" : "超限→需升级 L2") << std::endl;

        std::ofstream f(out_dir + "/result3_线性化适用性.csv", std::ios::binary);
        put_bom(f);
        f << "指标,观测值,阈值,判定\n";
        f << "max|θf|(rad)," << std::setprecision(12) << max_tf << "," << ANGLE_LIMIT << ","
          << (max_tf <= ANGLE_LIMIT ? "PASS" : "FAIL") << "\n";
        f << "max|θo|(rad)," << max_to << "," << ANGLE_LIMIT << ","
          << (max_to <= ANGLE_LIMIT ? "PASS" : "FAIL") << "\n";
        f << "max|u|/ℓ_e," << slip_ratio << "," << SLIP_LIMIT << ","
          << (slip_ratio <= SLIP_LIMIT ? "PASS" : "FAIL") << "\n";
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_线性化适用性.csv" << std::endl;
    }

    // ======== 响应幅值摘要 ========
    {
        double max_zf = 0, max_tf = 0, max_zo = 0, max_to = 0;
        double max_vf = 0, max_wf = 0, max_vo = 0, max_wo = 0;
        for (auto& s : fr.samples) {
            if (std::fabs(s.zf) > max_zf) max_zf = std::fabs(s.zf);
            if (std::fabs(s.vf) > max_vf) max_vf = std::fabs(s.vf);
            if (std::fabs(s.tf) > max_tf) max_tf = std::fabs(s.tf);
            if (std::fabs(s.wf) > max_wf) max_wf = std::fabs(s.wf);
            if (std::fabs(s.zo) > max_zo) max_zo = std::fabs(s.zo);
            if (std::fabs(s.vo) > max_vo) max_vo = std::fabs(s.vo);
            if (std::fabs(s.to) > max_to) max_to = std::fabs(s.to);
            if (std::fabs(s.wo) > max_wo) max_wo = std::fabs(s.wo);
        }
        std::ofstream f(out_dir + "/result3_响应幅值摘要.csv", std::ios::binary);
        put_bom(f);
        f << "运动量,中文含义,最大绝对值,单位\n";
        f << "zf,浮子垂荡位移," << std::setprecision(10) << max_zf << ",m\n";
        f << "dzf/dt,浮子垂荡速度," << max_vf << ",m/s\n";
        f << "θf,浮子纵摇角位移," << max_tf << ",rad\n";
        f << "dθf/dt,浮子纵摇角速度," << max_wf << ",rad/s\n";
        f << "zo,振子垂荡位移," << max_zo << ",m\n";
        f << "dzo/dt,振子垂荡速度," << max_vo << ",m/s\n";
        f << "θo,振子纵摇角位移," << max_to << ",rad\n";
        f << "dθo/dt,振子纵摇角速度," << max_wo << ",rad/s\n";
        f.close();
        std::cout << "\n====== 响应幅值摘要 ======" << std::endl;
        std::cout << "max|zf|=" << std::setprecision(6) << max_zf << " m  max|żf|=" << max_vf << " m/s" << std::endl;
        std::cout << "max|θf|=" << max_tf << " rad  max|θ̇f|=" << max_wf << " rad/s" << std::endl;
        std::cout << "max|zo|=" << max_zo << " m  max|żo|=" << max_vo << " m/s" << std::endl;
        std::cout << "max|θo|=" << max_to << " rad  max|θ̇o|=" << max_wo << " rad/s" << std::endl;
        std::cout << "已写入 " << out_dir << "/result3_响应幅值摘要.csv" << std::endl;
    }

    std::cout << "\n======== Q3 v3 求解完成 ========" << std::endl;
    std::cout << "全部结果输出至: " << out_dir << std::endl;

    return 0;
}
