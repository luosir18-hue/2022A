// 第二问求解器 (C++ 详细注释版)。Owner: 罗懿 (L, 建模手).
// =============================== 问题描述 ===============================
// 问题二在Q1垂荡动力学基础上,将PTO阻尼参数作为决策变量,
// 在题目给定的闭区间内最大化PTO的稳态平均机械输出功率.
// 情形一: 常量阻尼 c ∈ [0, 100000] N·s/m 的一维有界优化.
// 情形二: 幂律阻尼 D = c0·|v_r|^p·v_r, c0∈[0,100000], p∈[0,1] 的二维优化.
//
// =============================== 数据来源 ===============================
// 附件3 (Sheet1): 匹配 ω=2.2143 s⁻¹ 行 → C列=附加质量,E列=兴波阻尼,G列=激励力
// 附件4 (Sheet1!A:B): 浮子/振子质量,半径,密度,重力,弹簧刚度
// Kh = ρgπR_f² 由水线面面积公式计算.
//
// =============================== 模型依据 ===============================
// 模型/Q2/v1/数学模型.md 式(Q2-1)~(Q2-8).
//   控制方程(Q2-1): 二自由度垂荡运动,浮子受波浪激励+浮子/振子弹簧-PTO耦合
//   瞬时功率(Q2-2): P(t)=D(v_r)·v_r (阻尼力×相对速度,恒≥0)
//   平均功率(Q2-3): 稳态长期时间平均,剔除初始瞬态
//   复幅值模型(Q2-4,Q2-5): 仅对常量阻尼,用于解析校验
//   优化模型(Q2-6,Q2-7): 在有界闭区间内最大化稳态平均功率
//
// =============================== 算法结构 ===============================
// [模块A: XLSX解析] 手工解析ZIP+sheet1.xml,提取物理参数(从附件文件).
//   └── load_params: 附件3匹配ω行→m_a,b_h,F; 附件4参数名→值; Kh公式计算.
//
// [模块B: 时域仿真] 四阶龙格-库塔法(RK4)求解二自由度ODE,积分步长0.001s.
//   ├── rhs: 右端函数,实现式(Q2-1)浮子/振子加速度.
//   ├── rk4_step: RK4单步积分,局部O(dt⁵)全局O(dt⁴).
//   ├── steady_avg_power: 预热n_warmup周期+梯形积分n_eval周期→平均功率.
//   ├── simulate_traj: 完整轨迹记录,采样间隔0.2s.
//   └── check_steady_convergence: 逐周期输出功率+状态值,双判据(功率/状态)诊断收敛.
//
// [模块C: 频域解析解] 常量阻尼的闭式解.
//   └── analytical_power: P_lin(c)=½cω²|ẑ_o−ẑ_f|², 由复幅值矩阵方程导出.
//
// [模块D: 黄金分割搜索] 一维单峰函数优化.
//   └── golden_section<Func>: 区间[a,b]→每次缩至φ≈0.618倍,至b-a≤tol止.
//
// [模块E: 情形一优化] 常量阻尼.
//   ├── 解析: golden_section(f_anal, 0, 100000, tol=1e-6)
//   ├── 时域: golden_section(f_td, 0, 100000, tol=5) 用作验证
//   ├── 输出: 结果摘要CSV, 功率曲线CSV(40点解析+时域双列),
//   │         最优轨迹CSV(60周期×0.2s), 稳态验证CSV(40周期逐周期功率+状态)
//
// [模块F: 情形二优化] 幂律阻尼 (粗网格+局部精细+最终评估).
//   ├── 粗网格: 16×11点, c0按(i/15)^2.5非线性分布(小阻尼区密集), 20暖+10评
//   ├── 局部精细: 粗网最优±25%c0/±0.12p, 9×9点, 25暖+15评
//   ├── 最终评估: 30暖+30评
//   ├── 输出: 结果摘要CSV(含边界最优检测),
//   │         粗网格曲面CSV(176点), 精细网格邻域CSV(81点),
//   │         最优轨迹CSV, 稳态验证CSV
//
// =============================== 编译运行 ===============================
// 编译: g++ -std=c++17 -O2 -o solve solve.cpp -lz
// 运行: ./solve [--data <数据目录>] [--out <输出目录>]
//   默认数据: 自动搜索 题目信息/A题 (适配多级相对路径)
//   默认输出: 结果（包括各种过程数据与审查）/运行输出/run001 (自动递增)
//   CSV编码: UTF-8 BOM (Excel/VSCode不乱码)
// =======================================================================


#include <cmath>
#include <complex>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

namespace {

const double PI = 3.14159265358979323846;
const double OMEGA_TARGET = 2.2143;
const double CI_LOW = 0.0, CI_HIGH = 100000.0;
const double P_LOW = 0.0, P_HIGH = 1.0;
const double GOLDEN_R = (std::sqrt(5.0) - 1.0) / 2.0;

struct Params {
    double m_f, m_o, m_a, b_h, F, rho, g, R_f, k, Kh;
};

struct Damp {
    int type;
    double c, c0, pp;
};

double damp_force(double vr, const Damp& d) {
    if (d.type == 1) return d.c * vr;
    return d.c0 * std::pow(std::fabs(vr), d.pp) * vr;
}

// ======================================================================
// XLSX reader — same as Q1, adapted for ω=2.2143
// ======================================================================
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
    p.R_f = need4("浮子底半径");
    p.m_o = need4("振子质量");
    p.rho = need4("海水的密度");
    p.g = need4("重力加速度");
    p.k = need4("弹簧刚度");
    p.Kh = p.rho * p.g * PI * p.R_f * p.R_f;

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
                p.m_a = get3("C");
                p.b_h = get3("E");
                p.F = get3("G");
                found = true;
                break;
            }
        }
    }
    if (!found) throw std::runtime_error("附件3 未找到 ω=2.2143 对应的 Q2 行");
    return p;
}

// ======================================================================
// Time-domain simulation
// ======================================================================
void rhs(double t, const double x[4], double dx[4], const Params& par, const Damp& damp) {
    double xi = x[2] - x[0];
    double vr = x[3] - x[1];
    double D = damp_force(vr, damp);
    double exc = par.F * std::cos(OMEGA_TARGET * t);
    dx[0] = x[1];
    dx[1] = (exc - par.b_h * x[1] - par.Kh * x[0] + par.k * xi + D) / (par.m_f + par.m_a);
    dx[2] = x[3];
    dx[3] = (-par.k * xi - D) / par.m_o;
}

void rk4_step(double& t, double x[4], double dt, const Params& par, const Damp& damp) {
    double k1[4], k2[4], k3[4], k4[4], xt[4];
    rhs(t, x, k1, par, damp);
    for (int i = 0; i < 4; ++i) xt[i] = x[i] + 0.5 * dt * k1[i];
    rhs(t + 0.5 * dt, xt, k2, par, damp);
    for (int i = 0; i < 4; ++i) xt[i] = x[i] + 0.5 * dt * k2[i];
    rhs(t + 0.5 * dt, xt, k3, par, damp);
    for (int i = 0; i < 4; ++i) xt[i] = x[i] + dt * k3[i];
    rhs(t + dt, xt, k4, par, damp);
    for (int i = 0; i < 4; ++i) x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    t += dt;
}

double steady_avg_power(const Params& par, const Damp& damp,
                        int n_warmup, int n_eval, double dt = 0.001) {
    double T = 2.0 * PI / OMEGA_TARGET;
    double t_warmup = n_warmup * T;
    double t_eval = n_eval * T;
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    double t = 0.0;

    while (t < t_warmup - dt * 0.5) rk4_step(t, x, dt, par, damp);

    double vr = x[3] - x[1];
    double prev_P = damp_force(vr, damp) * vr;
    double power_integral = 0.0;
    while (t < t_warmup + t_eval - dt * 0.5) {
        rk4_step(t, x, dt, par, damp);
        vr = x[3] - x[1];
        double P = damp_force(vr, damp) * vr;
        power_integral += 0.5 * (P + prev_P) * dt;
        prev_P = P;
    }
    return power_integral / t_eval;
}

struct Traj {
    std::vector<double> t, zf, vf, zo, vo, power;
};
Traj simulate_traj(const Params& par, const Damp& damp,
                   double t_end, double sample_dt = 0.2, double rk_dt = 0.001) {
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    double t = 0.0;
    long sample_interval = (long)std::lround(sample_dt / rk_dt);
    Traj tr;
    tr.t.push_back(0.0);
    tr.zf.push_back(0.0); tr.vf.push_back(0.0);
    tr.zo.push_back(0.0); tr.vo.push_back(0.0);
    double vr0 = 0.0;
    tr.power.push_back(damp_force(vr0, damp) * vr0);
    long step = 0;
    while (t < t_end - rk_dt * 0.5) {
        rk4_step(t, x, rk_dt, par, damp);
        ++step;
        if (step % sample_interval == 0) {
            double ts = step / sample_interval * sample_dt;
            double vr = x[3] - x[1];
            tr.t.push_back(ts);
            tr.zf.push_back(x[0]); tr.vf.push_back(x[1]);
            tr.zo.push_back(x[2]); tr.vo.push_back(x[3]);
            tr.power.push_back(damp_force(vr, damp) * vr);
        }
    }
    return tr;
}

void check_steady_convergence(const Params& par, const Damp& damp,
                              int n_total, std::ostream& log) {
    double T = 2.0 * PI / OMEGA_TARGET;
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    double t = 0.0, dt = 0.001;
    double prev_cycle_power = 0.0;
    double prev_cycle_x0 = 0.0, prev_cycle_x2 = 0.0;
    int converged_power = -1, converged_state = -1;
    log << "周期,平均功率(W),周期末浮子位移(m),周期末振子位移(m)\n";
    for (int cyc = 1; cyc <= n_total; ++cyc) {
        double cyc_start = t;
        double cyc_integral = 0.0;
        double vr = x[3] - x[1];
        double prev_P = damp_force(vr, damp) * vr;
        while (t < cyc_start + T - dt * 0.5) {
            rk4_step(t, x, dt, par, damp);
            vr = x[3] - x[1];
            double P = damp_force(vr, damp) * vr;
            cyc_integral += 0.5 * (P + prev_P) * dt;
            prev_P = P;
        }
        double cyc_avg = cyc_integral / T;
        log << cyc << "," << std::setprecision(10) << cyc_avg << ","
            << x[0] << "," << x[2] << "\n";
        if (cyc > 1) {
            double rel_power = std::fabs(cyc_avg - prev_cycle_power)
                              / std::max(std::fabs(cyc_avg), 1e-12);
            if (rel_power < 1e-4 && converged_power < 0) converged_power = cyc;
            double rel_state = std::max(
                std::fabs(x[0] - prev_cycle_x0) / std::max(std::fabs(x[0]), 1e-8),
                std::fabs(x[2] - prev_cycle_x2) / std::max(std::fabs(x[2]), 1e-8));
            if (rel_state < 1e-3 && converged_state < 0) converged_state = cyc;
        }
        prev_cycle_power = cyc_avg;
        prev_cycle_x0 = x[0];
        prev_cycle_x2 = x[2];
    }
    if (converged_power > 0)
        log << "稳态成立周期_功率判据(|ΔP|/P<1e-4)," << converged_power << "\n";
    else
        log << "稳态成立周期_功率判据,未在" << n_total << "周期内收敛\n";
    if (converged_state > 0)
        log << "稳态成立周期_状态判据(|Δz|/|z|<1e-3)," << converged_state << "\n";
    else
        log << "稳态成立周期_状态判据,未在" << n_total << "周期内收敛\n";
}

// ======================================================================
// Analytical solution for constant damping
// ======================================================================
double analytical_power(double c, const Params& par) {
    using cd = std::complex<double>;
    double om2 = OMEGA_TARGET * OMEGA_TARGET;
    double kh = par.Kh, k = par.k, mf = par.m_f, mo = par.m_o;
    double ma = par.m_a, bh = par.b_h;

    cd Z11(-om2 * (mf + ma) + (kh + k), OMEGA_TARGET * (bh + c));
    cd Z12(-k, -OMEGA_TARGET * c);
    cd Z21(-k, -OMEGA_TARGET * c);
    cd Z22(-om2 * mo + k, OMEGA_TARGET * c);

    cd det = Z11 * Z22 - Z12 * Z21;
    cd zf = (Z22 * par.F) / det;
    cd zo = (-Z21 * par.F) / det;
    double amp_rel = std::abs(zo - zf);
    return 0.5 * c * om2 * amp_rel * amp_rel;
}

// ======================================================================
// 1D golden-section search
// ======================================================================
struct Opt1D {
    double x_opt, f_opt;
    std::vector<double> xs, fs;
};

template<typename Func>
Opt1D golden_section(Func f, double a, double b, double tol = 1e-6) {
    Opt1D res;
    double c = b - GOLDEN_R * (b - a);
    double d = a + GOLDEN_R * (b - a);
    double fc = f(c), fd = f(d);
    res.xs = {a, c, d, b};
    res.fs = {f(a), fc, fd, f(b)};
    int iter = 0;
    while (b - a > tol && iter < 200) {
        if (fc > fd) {
            b = d; d = c; fd = fc;
            c = b - GOLDEN_R * (b - a);
            fc = f(c);
        } else {
            a = c; c = d; fc = fd;
            d = a + GOLDEN_R * (b - a);
            fd = f(d);
        }
        ++iter;
    }
    double xm = (a + b) / 2.0;
    res.x_opt = xm;
    res.f_opt = f(xm);
    return res;
}

// ======================================================================
// 2D grid search + local refinement (power-law damping)
// ======================================================================
struct Point2D {
    double c0, pp, fval;
};

// ======================================================================
// Output helpers
// ======================================================================
void put_bom(std::ofstream& out) { out << "\xEF\xBB\xBF"; }

// ======================================================================
// Path resolution
// ======================================================================
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

}  // namespace

// ======================================================================
// main
// ======================================================================
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

    Params par = load_params(data_dir + "/附件3.xlsx", data_dir + "/附件4.xlsx");
    std::cout << "参数加载（从附件文件读取）：" << std::endl;
    std::cout << "  m_f=" << par.m_f << " m_o=" << par.m_o << " R_f=" << par.R_f << std::endl;
    std::cout << "  ρ=" << par.rho << " g=" << par.g << " k=" << par.k << " Kh=" << par.Kh << std::endl;
    std::cout << "  ω=" << OMEGA_TARGET << " m_a=" << par.m_a << " b_h=" << par.b_h
              << " F=" << par.F << std::endl;

    // ======= 常量阻尼一维优化 =======
    std::cout << "\n=== 情形一：常量阻尼 c ∈ [0, 100000] ===" << std::endl;

    auto f_anal = [&](double c) { return analytical_power(c, par); };
    auto opt_anal = golden_section(f_anal, CI_LOW, CI_HIGH, 1e-6);
    std::cout << "解析最优: c*=" << std::setprecision(10) << opt_anal.x_opt
              << " N·s/m  P*=" << opt_anal.f_opt << " W" << std::endl;

    auto f_td = [&](double c) -> double {
        Damp d; d.type = 1; d.c = c;
        return steady_avg_power(par, d, 25, 15);
    };
    auto opt_td = golden_section(f_td, CI_LOW, CI_HIGH, 5.0);
    std::cout << "时域最优: c*=" << std::setprecision(10) << opt_td.x_opt
              << " N·s/m  P*=" << opt_td.f_opt << " W" << std::endl;

    // 功率曲线: 40 points
    {
        std::ofstream f(out_dir + "/result2-1_功率曲线.csv", std::ios::binary);
        put_bom(f);
        f << "阻尼系数c(N·s/m),平均功率_解析(W),平均功率_时域(W)\n";
        int npt = 40;
        for (int i = 0; i <= npt; ++i) {
            double c = CI_LOW + (CI_HIGH - CI_LOW) * i / npt;
            double pan = analytical_power(c, par);
            Damp d; d.type = 1; d.c = c;
            double ptd = steady_avg_power(par, d, 20, 15);
            f << std::setprecision(10) << c << "," << pan << "," << ptd << "\n";
        }
        f.close();
    }

    // 最优参数时域详细轨迹
    {
        double c_opt = opt_anal.x_opt;
        Damp d; d.type = 1; d.c = c_opt;
        double T = 2.0 * PI / OMEGA_TARGET;
        double t_end = 60.0 * T;
        auto tr = simulate_traj(par, d, t_end);
        std::ofstream f(out_dir + "/result2-1_最优轨迹.csv", std::ios::binary);
        put_bom(f);
        f << "时间(s),浮子位移(m),浮子速度(m/s),振子位移(m),振子速度(m/s),瞬时功率(W)\n";
        for (size_t i = 0; i < tr.t.size(); ++i) {
            f << std::setprecision(10) << tr.t[i] << "," << tr.zf[i] << "," << tr.vf[i] << ","
              << tr.zo[i] << "," << tr.vo[i] << "," << tr.power[i] << "\n";
        }
        f.close();

        // 稳态收敛检查
        std::ofstream flog(out_dir + "/result2-1_稳态验证.csv", std::ios::binary);
        put_bom(flog);
        check_steady_convergence(par, d, 40, flog);
        flog.close();
    }

    // 最优结果摘要
    {
        std::ofstream f(out_dir + "/result2-1_最优.csv", std::ios::binary);
        put_bom(f);
        f << "项目,值\n";
        f << "最优阻尼系数c_opt(N·s/m)," << std::setprecision(10) << opt_anal.x_opt << "\n";
        f << "最大平均功率_解析(W)," << opt_anal.f_opt << "\n";
        f << "最大平均功率_时域(W)," << opt_td.f_opt << "\n";
        f << "解析与时域功率差值(W),"
          << std::fabs(opt_anal.f_opt - opt_td.f_opt) << "\n";
        f << "优化方法_解析,黄金分割搜索(容差1e-6)\n";
        f << "优化方法_时域,黄金分割搜索(容差5N·s/m)\n";
        f.close();
    }

    std::cout << "情形一完成。输出至 " << out_dir << std::endl;

    // ======= 幂律阻尼二维优化 =======
    std::cout << "\n=== 情形二：幂律阻尼 c0 ∈ [0, 100000], p ∈ [0, 1] ===" << std::endl;

    int grid_c0 = 16, grid_p = 11;
    std::vector<Point2D> grid;
    for (int i = 0; i < grid_c0; ++i) {
        double c0;
        if (i == 0) c0 = 0.0;
        else c0 = CI_LOW + (CI_HIGH - CI_LOW) * std::pow((double)i / (grid_c0 - 1), 2.5);
        c0 = std::min(c0, CI_HIGH);
        for (int j = 0; j < grid_p; ++j) {
            double pp = P_LOW + (P_HIGH - P_LOW) * j / (grid_p - 1);
            Damp d; d.type = 2; d.c0 = c0; d.pp = pp;
            double pavg = steady_avg_power(par, d, 20, 10);
            grid.push_back({c0, pp, pavg});
        }
    }

    Point2D best{0, 0, -1};
    for (auto& g : grid)
        if (g.fval > best.fval) best = g;
    std::cout << "粗网格最优: c0=" << best.c0 << " p=" << best.pp
              << " P=" << best.fval << " W" << std::endl;

    // 局部精细搜索
    std::vector<Point2D> neighbors;
    double sr_c0 = std::max(best.c0 * 0.25, 500.0);
    double sr_p = 0.12;
    int fine_pts = 9;
    for (int i = 0; i < fine_pts; ++i) {
        double fc0 = best.c0 - sr_c0 + 2.0 * sr_c0 * i / (fine_pts - 1);
        fc0 = std::max(P_LOW, std::min(CI_HIGH, fc0));
        for (int j = 0; j < fine_pts; ++j) {
            double fpp = best.pp - sr_p + 2.0 * sr_p * j / (fine_pts - 1);
            fpp = std::max(P_LOW, std::min(P_HIGH, fpp));
            Damp d; d.type = 2; d.c0 = fc0; d.pp = fpp;
            double pavg = steady_avg_power(par, d, 25, 15);
            neighbors.push_back({fc0, fpp, pavg});
        }
    }
    Point2D refined = best;
    for (auto& n : neighbors)
        if (n.fval > refined.fval) refined = n;

    // 输出精细网格搜索的邻域证据
    {
        std::ofstream f(out_dir + "/result2-2_精细网格_邻域搜索.csv", std::ios::binary);
        put_bom(f);
        f << "阻尼比例系数c0,幂指数p(无量纲),平均功率P(W)\n";
        for (auto& n : neighbors) {
            f << std::setprecision(10) << n.c0 << "," << n.pp << ","
              << n.fval << "\n";
        }
        f.close();
    }

    // 最终精细评估
    {
        Damp df; df.type = 2; df.c0 = refined.c0; df.pp = refined.pp;
        double final_p = steady_avg_power(par, df, 30, 30);
        refined.fval = final_p;
    }
    std::cout << "精细最优: c0*=" << std::setprecision(10) << refined.c0
              << " p*=" << refined.pp << " P*=" << refined.fval << " W" << std::endl;

    // 最优轨迹
    {
        Damp d; d.type = 2; d.c0 = refined.c0; d.pp = refined.pp;
        double T = 2.0 * PI / OMEGA_TARGET;
        double t_end = 60.0 * T;
        auto tr = simulate_traj(par, d, t_end);
        std::ofstream f(out_dir + "/result2-2_最优轨迹.csv", std::ios::binary);
        put_bom(f);
        f << "时间(s),浮子位移(m),浮子速度(m/s),振子位移(m),振子速度(m/s),瞬时功率(W)\n";
        for (size_t i = 0; i < tr.t.size(); ++i) {
            f << std::setprecision(10) << tr.t[i] << "," << tr.zf[i] << "," << tr.vf[i] << ","
              << tr.zo[i] << "," << tr.vo[i] << "," << tr.power[i] << "\n";
        }
        f.close();
    }

    // 稳态收敛
    {
        Damp d; d.type = 2; d.c0 = refined.c0; d.pp = refined.pp;
        std::ofstream flog(out_dir + "/result2-2_稳态验证.csv", std::ios::binary);
        put_bom(flog);
        check_steady_convergence(par, d, 40, flog);
        flog.close();
    }

    // 网格结果输出
    {
        std::ofstream f(out_dir + "/result2-2_功率曲面网格.csv", std::ios::binary);
        put_bom(f);
        f << "阻尼比例系数c0,幂指数p(无量纲),平均功率P(W)\n";
        for (auto& g : grid) {
            f << std::setprecision(10) << g.c0 << "," << g.pp << "," << g.fval << "\n";
        }
        f.close();
    }

    // 最优摘要
    {
        std::ofstream f(out_dir + "/result2-2_最优.csv", std::ios::binary);
        put_bom(f);
        f << "项目,值\n";
        f << "最优比例系数c0*," << std::setprecision(10) << refined.c0 << "\n";
        f << "最优幂指数p*," << refined.pp << "\n";
        f << "最大平均功率P*(W)," << refined.fval << "\n";
        f << "优化方法,粗网格" << grid_c0 << "×" << grid_p
          << "+局部精细" << fine_pts << "×" << fine_pts << "搜索\n";
        f << "粗网格最优c0," << best.c0 << "\n";
        f << "粗网格最优p," << best.pp << "\n";
        f << "粗网格最优P(W)," << best.fval << "\n";
        f << "最终评估周期_预热," << 30 << "\n";
        f << "最终评估周期_平均," << 30 << "\n";
        f << "是否边界最优,";
        if (refined.c0 >= CI_HIGH - 1e-6) f << "是(c0*触及上界100000)\n";
        else if (refined.c0 <= P_LOW + 1e-6) f << "是(c0*触及下界0)\n";
        else if (refined.pp >= P_HIGH - 1e-9) f << "是(p*触及上界1)\n";
        else if (refined.pp <= P_LOW + 1e-9) f << "是(p*触及下界0)\n";
        else f << "否(内部最优)\n";
        f.close();
    }

    std::cout << "情形二完成。输出至 " << out_dir << std::endl;
    std::cout << "\n全部 Q2 求解完成。" << std::endl;

    return 0;
}
