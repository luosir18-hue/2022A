// Q3 v1 求解器（C++ 详细中文注释版）. Owner: 罗懿 (L, 建模手).
// 数据来源: 题目信息/A题/附件3.xlsx (ω=1.7152 对应行), 附件4.xlsx (全部数值从文件读取, 无手写数据)
// 模型依据: 模型/Q3/v2/数学模型.md 式(Q3-1)~(Q3-24)
//
// ======================== 系统概述 ========================
// 本程序求解问题三的四自由度垂荡—纵摇耦合非线性动力系统初值问题。
// 采用拉格朗日力学推导的质量矩阵、速度二次项、势能梯度和耗散力，
// 以 4 阶显式 Runge-Kutta 法 (RK4) 在时域积分 40 个波浪周期的运动响应。
//
// 广义坐标: q = (zf, θf, ξ, θo)ᵀ
//   zf  — 浮子质心相对静平衡的垂荡位移 (m, 向上为正)
//   θf  — 浮子纵摇角 (rad)
//   ξ   — 振子沿中轴相对静平衡的滑移 (m)
//   θo  — 中轴与振子的绝对纵摇角 (rad)
// 一阶状态: y = (qᵀ, q̇ᵀ)ᵀ = (zf, θf, ξ, θo, żf, θ̇f, ξ̇, θ̇o)ᵀ 共 8 个分量
//
// 控制方程 (式 Q3-10 ~ Q3-16):
//   M(q)·q̈ + h(q,q̇) + ∇V(q) + ∂R/∂q̇ = Q(t)
//   q̈ = M(q)⁻¹ [Q(t) − h(q,q̇) − ∇V(q) − ∂R/∂q̇]
//   其中:
//   (Q3-11) 质量矩阵 M(q)    — 4×4 对称, 含位置依赖项
//   (Q3-12) 速度二次项 h      — 4 分量, 含 Christoffel 符号贡献
//   (Q3-13) 势能梯度 ∇V       — 4 分量
//   (Q3-14) 耗散广义力 ∂R/∂q̇ — 4 分量
//   Q(t) = (F·cosωt, L·cosωt, 0, 0)ᵀ  — 波浪激励
//
// 输出转换:
//   (Q3-22) 振子绝对垂荡位移 z_o^out = zf + d(cosθf−1) + ℓ·cosθo − ℓe
//   (Q3-23) 振子绝对垂荡速度 ż_o^out = żf − d·sinθf·θ̇f + ξ̇·cosθo − ℓ·sinθo·θ̇o
//
// 输出文件 (全部写入本版本运行输出目录):
//   result3.xlsx                — 正式结果 (9 列完整时序, 与 result3.xlsx 列序一致)
//   result3.csv                 — 同一数据 CSV 版
//   result3_key.csv             — 论文五时刻 (t=10,20,40,60,100 s)
//   result3_弹簧长度检查.csv     — 弹簧长度全程记录 (硬约束 ξ>-0.2019575)
//   result3_质量矩阵正定性.csv   — 质量矩阵最小特征值全程记录 (硬约束 M 正定)
//   result3_能量恒等式.csv       — 能量守恒对账 (式 Q3-18)
//   result3_线性模型对照.csv     — 线性模型对照 (诊断, 式 Q3-17)
//   result3_响应幅值摘要.csv     — 八个运动量最大幅值摘要

// ========== 头文件引用 ==========
#include <cmath>        // sin, cos, atan2, fabs, sqrt 等数学函数
#include <cstdint>      // uint8_t, uint16_t, uint32_t 固定宽度整数
#include <cctype>       // isalpha 字符判断
#include <filesystem>   // 路径解析与目录创建
#include <fstream>      // 文件读写
#include <iomanip>      // setprecision, setw 输出格式化
#include <iostream>     // cout, cerr 控制台输出
#include <map>          // 键值对 (XLSX 单元格存储)
#include <sstream>      // ostringstream 字符串流
#include <stdexcept>    // runtime_error 异常
#include <string>       // std::string
#include <vector>       // 动态数组
#include <zlib.h>       // zlib inflate/deflate (XLSX 读写压缩)

namespace {  // 匿名命名空间: 内部符号限定在本文件, 避免跨文件 ODR 冲突

// ========== 全局常数 ==========
const double PI = 3.14159265358979323846;  // 圆周率 (高精度)
const double OMEGA_TARGET = 1.7152;        // 题面给定目标圆频率, 用于在附件3中定位对应行 (s⁻¹)
const double C_L = 10000.0;                // 题面给定直线 PTO 阻尼系数 (N·s/m)
const double C_THETA = 1000.0;             // 题面给定旋转 PTO 阻尼系数 (N·m·s/rad)

// ========== 参数结构体 ==========
// 存储从附件3/附件4读取的全部参数及派生锁定常数
struct Params {
    // -- 附件4: 结构/几何/弹簧参数 --
    double m_f, m_o;            // 浮子质量、振子质量 (kg)
    double R_f, H_c, H_n;       // 浮子底半径、圆柱高、圆锥高 (m)
    double R_o, H_o;            // 振子半径、高度 (m)
    double rho, g;              // 海水密度 (kg/m³)、重力加速度 (m/s²)
    double k, l0;               // 直线弹簧刚度 (N/m)、原长 (m)
    double k_theta, K_theta;    // 扭转弹簧刚度 (N·m/rad)、静水恢复力矩系数 (N·m/rad)

    // -- 附件3: 水动力系数 (ω=1.7152 对应行) --
    double omega;               // 实际圆频率 (从附件3 B 列读取)
    double m_a;                 // 垂荡附加质量 (kg)
    double J_a;                 // 纵摇附加转动惯量 (kg·m²)
    double b_h;                 // 垂荡兴波阻尼 (N·s/m)
    double b_theta;             // 纵摇兴波阻尼 (N·m·s)
    double F_amp;               // 波浪激励力振幅 (N)
    double L_amp;               // 波浪激励力矩振幅 (N·m)

    // -- 派生锁定常数 (load_params 一次算得后固定, 不再重算) --
    double Kh;        // 垂荡静水恢复刚度 (N/m), Kh = ρ·g·π·R_f²
    double z_G;       // 浮子壳体质心高度 (m)
    double d;         // 转轴相对浮子质心的轴向偏置 (m), d<0 表示转轴在质心下方
    double I_f;       // 浮子绕质心横轴的结构转动惯量 (kg·m²)
    double I_o;       // 振子绕自身质心横轴的转动惯量 (kg·m²)
    double l_e;       // 弹簧静平衡长度 (m), 振子质心到转轴的平衡距离
};

// ====================================================================
// 【模块 A】XLSX 读取器
// 功能: 手写 ZIP 解析 + XML 提取, 只依赖 zlib, 无需任何 Excel 库。
//       用于读取附件3.xlsx 和附件4.xlsx 中的数值。
// ====================================================================

// 以二进制读取整个文件到字符串
std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("无法打开文件: " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// zlib 原始 deflate 解压 (ZIP 的 DEFLATE 压缩条目)
// MAX_WBITS 取负值表示不解析 zlib 头, 直接处理原始 deflate 流
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

// 单元格: 字符串或数值
struct Cell { bool is_str = false; double num = 0.0; std::string str; };

// 解析 xl/sharedStrings.xml 得到共享字符串表
// XLSX 将字符串按索引存储, 单元格通过索引引用
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

// 解析 xl/worksheets/sheet1.xml, 返回单元格引用 (如 "A1") → Cell 的映射
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

// ZIP 解析: 从文件尾部 EOCD 定位中央目录, 遍历条目解压得到各 XML 文件
std::map<std::string, std::string> unzip_entries(const std::string& data) {
    std::map<std::string, std::string> out;
    // 反向搜索 EOCD 签名 0x50 0x4b 0x05 0x06 ("PK\x05\x06")
    size_t eocd = std::string::npos;
    for (size_t i = data.size() >= 22 ? data.size() - 22 : 0; ; --i) {
        if (i + 4 <= data.size() &&
            (uint8_t)data[i] == 0x50 && (uint8_t)data[i + 1] == 0x4b &&
            (uint8_t)data[i + 2] == 0x05 && (uint8_t)data[i + 3] == 0x06) { eocd = i; break; }
        if (i == 0) break;
    }
    if (eocd == std::string::npos) throw std::runtime_error("不是有效的 xlsx(zip) 文件");
    // EOCD 偏移量 16 处的 4 字节 = 中央目录起始偏移
    uint32_t co = (uint32_t)(uint8_t)data[eocd + 16]
        | ((uint32_t)(uint8_t)data[eocd + 17] << 8)
        | ((uint32_t)(uint8_t)data[eocd + 18] << 16)
        | ((uint32_t)(uint8_t)data[eocd + 19] << 24);
    size_t p = co;
    // 遍历中央目录条目 (签名 0x50 0x4b 0x01 0x02)
    while (p + 4 <= data.size() &&
           (uint8_t)data[p] == 0x50 && (uint8_t)data[p + 1] == 0x4b &&
           (uint8_t)data[p + 2] == 0x01 && (uint8_t)data[p + 3] == 0x02) {
        uint16_t m = (uint16_t)((uint8_t)data[p + 10] | ((uint8_t)data[p + 11] << 8));  // 压缩方法
        uint32_t cs = (uint32_t)(uint8_t)data[p + 20]                                  // 压缩后大小
            | ((uint32_t)(uint8_t)data[p + 21] << 8)
            | ((uint32_t)(uint8_t)data[p + 22] << 16)
            | ((uint32_t)(uint8_t)data[p + 23] << 24);
        uint16_t fnl = (uint16_t)((uint8_t)data[p + 28] | ((uint8_t)data[p + 29] << 8)); // 文件名长度
        uint16_t xl = (uint16_t)((uint8_t)data[p + 30] | ((uint8_t)data[p + 31] << 8));  // 扩展字段长度
        uint16_t cl = (uint16_t)((uint8_t)data[p + 32] | ((uint8_t)data[p + 33] << 8));  // 注释长度
        uint32_t lo = (uint32_t)(uint8_t)data[p + 42]                                   // 本地文件头偏移
            | ((uint32_t)(uint8_t)data[p + 43] << 8)
            | ((uint32_t)(uint8_t)data[p + 44] << 16)
            | ((uint32_t)(uint8_t)data[p + 45] << 24);
        std::string name = data.substr(p + 46, fnl);
        size_t lp = lo;
        uint16_t lfnl = (uint16_t)((uint8_t)data[lp + 26] | ((uint8_t)data[lp + 27] << 8));
        uint16_t lxl = (uint16_t)((uint8_t)data[lp + 28] | ((uint8_t)data[lp + 29] << 8));
        size_t ds = lp + 30 + lfnl + lxl;  // 压缩数据起始
        std::string comp = data.substr(ds, cs);
        out[name] = (m == 0) ? comp : inflate_raw(comp);  // 方法0=存储, 否则解压
        p += 46 + fnl + xl + cl;
    }
    return out;
}

// 读取整个 xlsx 文件, 返回 sheet1 的单元格映射
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

// 从 "A1" 格式引用提取列字母
std::string col_letters(const std::string& ref) {
    size_t i = 0;
    while (i < ref.size() && std::isalpha((unsigned char)ref[i])) ++i;
    return ref.substr(0, i);
}
// 从 "A1" 格式引用提取行号
int row_number(const std::string& ref) {
    size_t i = 0;
    while (i < ref.size() && std::isalpha((unsigned char)ref[i])) ++i;
    return std::stoi(ref.substr(i));
}
// 去掉参数名末尾的 " (单位)" 后缀, 如 "浮子质量 (kg)" → "浮子质量"
std::string norm_label(const std::string& s) {
    size_t pos = s.rfind(" (");
    if (pos != std::string::npos) return s.substr(0, pos);
    return s;
}

// ====================================================================
// 【模块 B】参数加载 load_params
// 功能: 从附件4读取结构参数, 从附件3读取 ω=1.7152 行的水动力系数,
//       并一次性计算派生锁定常数 (Kh, z_G, d, I_f, I_o, l_e)。
//       所有数值均来自文件, 不在代码中手写任何数据。
// ====================================================================
Params load_params(const std::string& attach3, const std::string& attach4) {
    // 附件4: A 列为参数名(含单位), B 列为数值
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

    // 附件3: B 列为频率 ω, 定位 ω=1.7152 行, 读取实际频率及六个水动力系数
    auto c3 = read_xlsx(attach3);
    bool found = false;
    for (auto& kv : c3) {
        if (col_letters(kv.first) == "B" && !kv.second.is_str) {
            int r = row_number(kv.first);
            if (std::fabs(kv.second.num - OMEGA_TARGET) < 5e-4) {  // 频率容差匹配
                auto get3 = [&](const std::string& col) -> double {
                    std::string ref = col + std::to_string(r);
                    auto it = c3.find(ref);
                    if (it == c3.end() || it->second.is_str)
                        throw std::runtime_error("附件3 缺少列 " + col);
                    return it->second.num;
                };
                p.omega = get3("B");   // 实际圆频率 (从文件读取)
                p.m_a = get3("C");     // 垂荡附加质量
                p.J_a = get3("D");     // 纵摇附加转动惯量
                p.b_h = get3("E");     // 垂荡兴波阻尼
                p.b_theta = get3("F"); // 纵摇兴波阻尼
                p.F_amp = get3("G");   // 波浪激励力振幅
                p.L_amp = get3("H");   // 波浪激励力矩振幅
                found = true;
                break;
            }
        }
    }
    if (!found) throw std::runtime_error("附件3 未找到 ω=1.7152 对应的 Q3 行");

    // ---- 派生锁定常数 (式 Q3-1 至 Q3-5), 一次算得后固定, 不再重算 ----
    p.Kh = p.rho * p.g * PI * p.R_f * p.R_f;  // 垂荡静水恢复刚度 (Q3-1 前)

    double S_t = PI * p.R_f * p.R_f;          // 顶部圆面积
    double S_c = 2.0 * PI * p.R_f * p.H_c;    // 圆柱侧面积
    double S_n = PI * p.R_f * std::sqrt(p.R_f * p.R_f + p.H_n * p.H_n); // 圆锥侧面积
    double S_total = S_t + S_c + S_n;
    double sigma = p.m_f / S_total;           // 统一面密度 (kg/m²)
    double m_t = sigma * S_t;                 // 顶部圆面质量
    double m_c = sigma * S_c;                 // 圆柱侧壳质量
    double m_n = sigma * S_n;                 // 圆锥侧壳质量
    double z_t = p.H_n + p.H_c;               // 顶部质心高度 (圆锥尖端为轴向零点)
    double z_c = p.H_n + p.H_c * 0.5;         // 圆柱质心高度
    double z_n = (2.0 / 3.0) * p.H_n;         // 圆锥质心高度
    p.z_G = (m_t * z_t + m_c * z_c + m_n * z_n) / p.m_f;  // 浮子质心高度 (Q3-1)
    p.d = p.H_n - p.z_G;                      // 转轴偏置 (Q3-2)
    p.I_f = m_t * (p.R_f * p.R_f / 4.0 + (z_t - p.z_G) * (z_t - p.z_G))
           + m_c * (p.R_f * p.R_f / 2.0 + p.H_c * p.H_c / 12.0 + (z_c - p.z_G) * (z_c - p.z_G))
           + m_n * (p.R_f * p.R_f / 4.0 + p.H_n * p.H_n / 18.0 + (z_n - p.z_G) * (z_n - p.z_G)); // 浮子转动惯量 (Q3-3)
    p.I_o = p.m_o / 12.0 * (3.0 * p.R_o * p.R_o + p.H_o * p.H_o);  // 振子转动惯量 (Q3-4)
    p.l_e = p.l0 + p.H_o * 0.5 - p.m_o * p.g / p.k;  // 弹簧平衡长度 (Q3-5)

    return p;
}

// ====================================================================
// 【模块 C】4×4 线性方程组求解 (带列主元高斯消元)
// 功能: 求解 M·q̈ = F 得到加速度 q̈
// ====================================================================
bool solve4x4(const double A[4][4], const double b[4], double x[4]) {
    double M[4][4], rhs[4];
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = A[i][j]; rhs[i] = b[i]; }

    // 前向消元: 逐列消去, 每列选主元
    for (int col = 0; col < 4; ++col) {
        int best = col;   // 列主元: 找绝对值最大的行
        for (int row = col + 1; row < 4; ++row)
            if (std::fabs(M[row][col]) > std::fabs(M[best][col])) best = row;
        if (std::fabs(M[best][col]) < 1e-14) return false;  // 主元过小 → 奇异
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
    // 回代求解
    for (int i = 3; i >= 0; --i) {
        double sum = rhs[i];
        for (int j = i + 1; j < 4; ++j) sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }
    return true;
}

// ====================================================================
// 【模块 D】Jacobi 特征值算法 (4×4 实对称矩阵)
// 功能: 计算质量矩阵 M 的全部特征值, 用于检验 M 是否正定 (硬约束 min_eig>0)
// ====================================================================
void jacobi_eigenvalues(double A[4][4], double ev[4]) {
    for (int iter = 0; iter < 30; ++iter) {
        // 找最大非对角元位置 (p,q)
        int p = 0, q = 1;
        double max_val = std::fabs(A[0][1]);
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                if (std::fabs(A[i][j]) > max_val) { max_val = std::fabs(A[i][j]); p = i; q = j; }
        if (max_val < 1e-14) break;  // 已对角化
        // Givens 旋转角 θ = ½·atan2(2·Apq, Aqq−App)
        double theta = 0.5 * std::atan2(2.0 * A[p][q], A[q][q] - A[p][p]);
        double c = std::cos(theta), s = std::sin(theta);
        // 对矩阵施加正交相似变换, 消去 Apq 和 Aqp
        for (int i = 0; i < 4; ++i) {
            if (i != p && i != q) {
                double aip = A[i][p], aiq = A[i][q];
                A[i][p] = c * aip - s * aiq;
                A[p][i] = A[i][p];
                A[i][q] = s * aip + c * aiq;
                A[q][i] = A[i][q];
            }
        }
        double app = A[p][p], aqq = A[q][q], apq = A[p][q];
        A[p][p] = c * c * app + s * s * aqq - 2.0 * c * s * apq;
        A[q][q] = s * s * app + c * c * aqq + 2.0 * c * s * apq;
        A[p][q] = A[q][p] = 0.0;
    }
    for (int i = 0; i < 4; ++i) ev[i] = A[i][i];
}

// ====================================================================
// 【模块 E】物理方程 — 实现式(Q3-11)~式(Q3-16)
// 这是整个求解器的核心, 直接翻译数学模型, 不补任何物理项。
// ====================================================================

// 计算质量矩阵 M(q) — 式(Q3-11)
// 输入 q[4]=(zf,θf,ξ,θo), 输出 M[4][4] 行优先存储
// M 对称, 位置依赖项使用瞬时值:
//   ℓ = ℓ_e + ξ  (振子质心到转轴瞬时距离)
//   diff = θo − θf (相对转角)
void compute_mass_matrix(const double q[4], const Params& p, double M[4][4]) {
    double l = p.l_e + q[2];                       // 瞬时距离
    double sf = std::sin(q[1]), so = std::sin(q[3]);
    double co = std::cos(q[3]);
    double diff = q[3] - q[1];
    double sd = std::sin(diff), cd = std::cos(diff);

    M[0][0] = p.m_f + p.m_a + p.m_o;
    M[0][1] = -p.m_o * p.d * sf;
    M[0][2] = p.m_o * co;
    M[0][3] = -p.m_o * l * so;

    M[1][0] = M[0][1];                            // 对称
    M[1][1] = p.I_f + p.J_a + p.m_o * p.d * p.d;
    M[1][2] = p.m_o * p.d * sd;
    M[1][3] = p.m_o * p.d * l * cd;

    M[2][0] = M[0][2];                            // 对称
    M[2][1] = M[1][2];                            // 对称
    M[2][2] = p.m_o;
    M[2][3] = 0.0;

    M[3][0] = M[0][3];                            // 对称
    M[3][1] = M[1][3];                            // 对称
    M[3][2] = 0.0;                                // 对称
    M[3][3] = p.I_o + p.m_o * l * l;
}

// 计算速度二次项 h(q,q̇) — 式(Q3-12)
// 即 (d/dt)(∂T/∂q̇) − ∂T/∂q 中不含 q̈ 的部分, 含 Christoffel 符号贡献
void compute_h(const double q[4], const double qd[4], const Params& p, double h[4]) {
    double stf = qd[1], sxi = qd[2], sto = qd[3];  // θ̇f, ξ̇, θ̇o
    double l = p.l_e + q[2];
    double cf = std::cos(q[1]), co = std::cos(q[3]);
    double so = std::sin(q[3]);
    double diff = q[3] - q[1];
    double sd = std::sin(diff), cd = std::cos(diff);

    h[0] = p.m_o * (-p.d * stf * stf * cf - l * sto * sto * co - 2.0 * sto * sxi * so);
    h[1] = p.m_o * p.d * (l * sto * sto * (-sd) + 2.0 * sto * sxi * cd);
    h[2] = p.m_o * (-p.d * stf * stf * cd - l * sto * sto);
    h[3] = p.m_o * l * (-p.d * stf * stf * (-sd) + 2.0 * sto * sxi);
}

// 计算势能梯度 ∇V(q) — 式(Q3-13)
void compute_gradV(const double q[4], const Params& p, double g[4]) {
    double l = p.l_e + q[2];
    g[0] = p.Kh * q[0];                                     // ∂V/∂zf
    g[1] = p.K_theta * q[1] + p.k_theta * (q[1] - q[3]) - p.m_o * p.g * p.d * std::sin(q[1]); // ∂V/∂θf
    g[2] = p.k * q[2] + p.m_o * p.g * (std::cos(q[3]) - 1.0);                                  // ∂V/∂ξ
    g[3] = p.k_theta * (q[3] - q[1]) - p.m_o * p.g * l * std::sin(q[3]);                       // ∂V/∂θo
}

// 计算耗散广义力 ∂R/∂q̇ — 式(Q3-14)
void compute_dissipation(const double qd[4], const Params& p, double d[4]) {
    d[0] = p.b_h * qd[0];                                 // 垂荡兴波阻尼
    d[1] = p.b_theta * qd[1] + C_THETA * (qd[1] - qd[3]); // 纵摇兴波阻尼 + 旋转 PTO
    d[2] = C_L * qd[2];                                   // 直线 PTO
    d[3] = C_THETA * (qd[3] - qd[1]);                     // 旋转 PTO
}

// 右端函数 rhs(t,y)→dy/dt — 式(Q3-10)~式(Q3-16) 的一阶系统
// dy[0..3] = q̇ (位置导数), dy[4..7] = q̈ (加速度)
// 加速度由解线性方程组 M·q̈ = Q − h − ∇V − ∂R/∂q̇ 得到
void rhs(double t, const double y[8], double dy[8], const Params& p) {
    const double* q = y;         // 前 4 分量: 广义坐标
    const double* qd = y + 4;    // 后 4 分量: 广义速度

    double M[4][4];
    compute_mass_matrix(q, p, M);

    double h[4];
    compute_h(q, qd, p, h);

    double gV[4];
    compute_gradV(q, p, gV);

    double dR[4];
    compute_dissipation(qd, p, dR);

    // 外部激励 Q(t) = (F·cosωt, L·cosωt, 0, 0)ᵀ
    double ext[4] = {p.F_amp * std::cos(p.omega * t), p.L_amp * std::cos(p.omega * t), 0.0, 0.0};

    // 组装力向量 F = Q − h − ∇V − ∂R/∂q̇
    double F_vec[4];
    for (int i = 0; i < 4; ++i)
        F_vec[i] = ext[i] - h[i] - gV[i] - dR[i];

    // 解 M·q̈ = F_vec
    double qdd[4];
    if (!solve4x4(M, F_vec, qdd)) {
        std::ostringstream ss;
        ss << "质量矩阵奇异 at t=" << t << " q=("
           << q[0] << "," << q[1] << "," << q[2] << "," << q[3] << ")";
        throw std::runtime_error(ss.str());
    }

    for (int i = 0; i < 4; ++i) dy[i] = qd[i];
    for (int i = 0; i < 4; ++i) dy[4 + i] = qdd[i];
}

// ====================================================================
// 【模块 F】RK4 8 状态单步积分
// 功能: 内部步长 dt=0.001 s (远小于输出间隔 0.2 s, 保证精度)
// ====================================================================
void rk4_step(double& t, double y[8], double dt, const Params& p) {
    double k1[8], k2[8], k3[8], k4[8], yt[8];
    rhs(t, y, k1, p);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
    rhs(t + 0.5 * dt, yt, k2, p);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
    rhs(t + 0.5 * dt, yt, k3, p);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + dt * k3[i];
    rhs(t + dt, yt, k4, p);
    for (int i = 0; i < 8; ++i) y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    t += dt;
}

// ====================================================================
// 【模块 G】输出转换 — 从广义坐标还原振子绝对垂荡运动
// ====================================================================
// 振子绝对垂荡位移 — 式(Q3-22)
double z_o_out(const double q[4], const Params& p) {
    double l = p.l_e + q[2];
    return q[0] + p.d * (std::cos(q[1]) - 1.0) + l * std::cos(q[3]) - p.l_e;
}
// 振子绝对垂荡速度 — 式(Q3-23)
double zd_o_out(const double q[4], const double qd[4], const Params& p) {
    double l = p.l_e + q[2];
    return qd[0] - p.d * std::sin(q[1]) * qd[1] + qd[2] * std::cos(q[3]) - l * std::sin(q[3]) * qd[3];
}

// 机械能 E = T + V — 用于能量恒等式验证式(Q3-18)
// T = ½(m_f+m_a)żf² + ½(I_f+J_a)θ̇f² + ½·m_o(ẋo²+żo²) + ½·I_o·θ̇o²
// V = ½Kh·zf² + ½k·ξ² + ½K_θ·θf² + ½k_θ(θo−θf)² + m_o·g[d(cosθf−1)+ℓ(cosθo−1)]
double compute_energy(const double q[4], const double qd[4], const Params& p) {
    double l = p.l_e + q[2];
    double sf = std::sin(q[1]), cf = std::cos(q[1]);
    double so = std::sin(q[3]), co = std::cos(q[3]);

    // 振子质心水平速度 ẋo = d·cosθf·θ̇f + ℓ·cosθo·θ̇o + ξ̇·sinθo
    double xod = p.d * cf * qd[1] + l * co * qd[3] + qd[2] * so;
    // 振子质心垂向速度 żo = żf − d·sinθf·θ̇f + ξ̇·cosθo − ℓ·sinθo·θ̇o
    double zod = qd[0] - p.d * sf * qd[1] + qd[2] * co - l * so * qd[3];

    double T = 0.5 * (p.m_f + p.m_a) * qd[0] * qd[0]
             + 0.5 * (p.I_f + p.J_a) * qd[1] * qd[1]
             + 0.5 * p.m_o * (xod * xod + zod * zod)
             + 0.5 * p.I_o * qd[3] * qd[3];

    double V = 0.5 * p.Kh * q[0] * q[0]
             + 0.5 * p.k * q[2] * q[2]
             + 0.5 * p.K_theta * q[1] * q[1]
             + 0.5 * p.k_theta * (q[3] - q[1]) * (q[3] - q[1])
             + p.m_o * p.g * (p.d * (cf - 1.0) + l * (co - 1.0));

    return T + V;
}

// 能量变化率解析式 — 式(Q3-18) 右端
// Ė = F·cosωt·żf + L·cosωt·θ̇f − b_h·żf² − b_θ·θ̇f² − c_l·ξ̇² − c_θ(θ̇o−θ̇f)²
double energy_rate_analytic(const double qd[4], const Params& p, double t) {
    double input = p.F_amp * std::cos(p.omega * t) * qd[0]
                 + p.L_amp * std::cos(p.omega * t) * qd[1];
    double dissip = p.b_h * qd[0] * qd[0]
                  + p.b_theta * qd[1] * qd[1]
                  + C_L * qd[2] * qd[2]
                  + C_THETA * (qd[3] - qd[1]) * (qd[3] - qd[1]);
    return input - dissip;
}

// ====================================================================
// 【模块 H】线性模型 (式 Q3-17) — 仅诊断用, 不写入 result3.xlsx
// 功能: 在平衡点冻结几何的线性化系统, 用于对照有限转角几何耦合强度
// ====================================================================
void linear_rhs(double t, const double y[8], double dy[8], const Params& p) {
    double M0[4][4] = {   // 线性质量矩阵 (q=0 处冻结)
        {p.m_f + p.m_a + p.m_o, 0, p.m_o, 0},
        {0, p.I_f + p.J_a + p.m_o * p.d * p.d, 0, p.m_o * p.d * p.l_e},
        {p.m_o, 0, p.m_o, 0},
        {0, p.m_o * p.d * p.l_e, 0, p.I_o + p.m_o * p.l_e * p.l_e}
    };
    double C0[4][4] = {   // 线性阻尼矩阵
        {p.b_h, 0, 0, 0},
        {0, p.b_theta + C_THETA, 0, -C_THETA},
        {0, 0, C_L, 0},
        {0, -C_THETA, 0, C_THETA}
    };
    double K0[4][4] = {   // 线性刚度矩阵
        {p.Kh, 0, 0, 0},
        {0, p.K_theta + p.k_theta - p.m_o * p.g * p.d, 0, -p.k_theta},
        {0, 0, p.k, 0},
        {0, -p.k_theta, 0, p.k_theta - p.m_o * p.g * p.l_e}
    };

    const double* q = y;
    const double* qd = y + 4;
    double ext[4] = {p.F_amp * std::cos(p.omega * t), p.L_amp * std::cos(p.omega * t), 0.0, 0.0};
    double F_vec[4];
    for (int i = 0; i < 4; ++i) {
        double sumC = 0.0, sumK = 0.0;
        for (int j = 0; j < 4; ++j) { sumC += C0[i][j] * qd[j]; sumK += K0[i][j] * q[j]; }
        F_vec[i] = ext[i] - sumC - sumK;
    }

    double qdd[4];
    if (!solve4x4(M0, F_vec, qdd)) {
        for (int i = 0; i < 8; ++i) dy[i] = 0.0;
        return;
    }
    for (int i = 0; i < 4; ++i) dy[i] = qd[i];
    for (int i = 0; i < 4; ++i) dy[4 + i] = qdd[i];
}

// 线性模型 RK4 单步积分
void rk4_step_linear(double& t, double y[8], double dt, const Params& p) {
    double k1[8], k2[8], k3[8], k4[8], yt[8];
    linear_rhs(t, y, k1, p);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
    linear_rhs(t + 0.5 * dt, yt, k2, p);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
    linear_rhs(t + 0.5 * dt, yt, k3, p);
    for (int i = 0; i < 8; ++i) yt[i] = y[i] + dt * k3[i];
    linear_rhs(t + dt, yt, k4, p);
    for (int i = 0; i < 8; ++i) y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    t += dt;
}

// ====================================================================
// 【模块 I】采样点数据结构
// ====================================================================
struct Sample {
    double t;        // 时间 (s)
    double zf, vf;   // 浮子垂荡位移 (m)、速度 (m/s)
    double tf, wf;   // 浮子纵摇角 (rad)、角速度 (rad/s)
    double zo, vo;   // 振子绝对垂荡位移 (m)、速度 (m/s)
    double to, wo;   // 振子绝对纵摇角 (rad)、角速度 (rad/s)
};

struct FullRec {
    std::vector<Sample> samples;          // 输出采样点 (间隔 0.2 s)
    std::vector<double> energy;           // 机械能 E (J)
    std::vector<double> energy_balance;   // E(0)+∫Ė_ana dt, 能量守恒对账 (J)
    std::vector<double> spring_length;    // 弹簧瞬时长度 ℓ (m)
    std::vector<double> min_eigenvalue;   // 质量矩阵最小特征值
};

// ====================================================================
// 【模块 J】正式仿真 simulate_full — 非线性 L2 模型 (式 Q3-16)
// 功能: RK4 积分 40 个波浪周期, 每 0.2 s 采样一次, 同时累计能量对账
// ====================================================================
FullRec simulate_full(const Params& p, double rk_dt = 0.001) {
    double T = 2.0 * PI / p.omega;         // 波浪周期
    double t_max = 40.0 * T;               // 仿真时长 40 周期
    double sample_dt = 0.2;                // 输出采样间隔
    long sample_interval = (long)std::lround(sample_dt / rk_dt); // 200 内部步 = 1 输出步

    double y[8] = {0};   // 初始状态: 全部为零 (静水平衡, 式 Q3-20)
    double t = 0.0;

    FullRec fr;
    fr.samples.push_back({0,0,0,0,0,0,0,0,0});  // t=0 记录
    double E0 = compute_energy(y, y + 4, p);
    fr.energy.push_back(E0);
    fr.energy_balance.push_back(E0);
    fr.spring_length.push_back(p.l_e);

    // 初始质量矩阵最小特征值
    double M_tmp[4][4], ev[4];
    compute_mass_matrix(y, p, M_tmp);
    jacobi_eigenvalues(M_tmp, ev);
    double min_ev = ev[0];
    for (int i = 1; i < 4; ++i) if (ev[i] < min_ev) min_ev = ev[i];
    fr.min_eigenvalue.push_back(min_ev);

    // 能量对账: 按内部步长对 Ė_ana 做梯形积分, 验证式(Q3-18)
    double balance = E0;
    double Edot_prev = energy_rate_analytic(y + 4, p, 0.0);

    long step = 0;
    while (t < t_max - rk_dt * 0.5) {
        rk4_step(t, y, rk_dt, p);
        ++step;

        // 累计能量对账
        double Edot_curr = energy_rate_analytic(y + 4, p, t);
        balance += 0.5 * (Edot_prev + Edot_curr) * rk_dt;
        Edot_prev = Edot_curr;

        // 每到输出时刻采样并记录诊断数据
        if (step % sample_interval == 0) {
            double ts = step / (double)sample_interval * sample_dt;
            Sample s;
            s.t = ts;
            s.zf = y[0]; s.vf = y[4];                 // 浮子垂荡
            s.tf = y[1]; s.wf = y[5];                 // 浮子纵摇
            s.zo = z_o_out(y, p);                     // 振子绝对垂荡 (式 Q3-22)
            s.vo = zd_o_out(y, y + 4, p);             // 振子绝对垂荡速度 (式 Q3-23)
            s.to = y[3]; s.wo = y[7];                 // 振子绝对纵摇
            fr.samples.push_back(s);

            fr.energy.push_back(compute_energy(y, y + 4, p));
            fr.energy_balance.push_back(balance);
            fr.spring_length.push_back(p.l_e + y[2]);  // ℓ = ℓe + ξ

            compute_mass_matrix(y, p, M_tmp);
            jacobi_eigenvalues(M_tmp, ev);
            min_ev = ev[0];
            for (int i = 1; i < 4; ++i) if (ev[i] < min_ev) min_ev = ev[i];
            fr.min_eigenvalue.push_back(min_ev);
        }
    }
    return fr;
}

// ====================================================================
// 【模块 K】线性模型仿真 simulate_linear (对照用)
// ====================================================================
std::vector<Sample> simulate_linear(const Params& p, double rk_dt = 0.001) {
    double T = 2.0 * PI / p.omega;
    double t_max = 40.0 * T;
    double sample_dt = 0.2;
    long sample_interval = (long)std::lround(sample_dt / rk_dt);

    double y[8] = {0};
    double t = 0.0;

    std::vector<Sample> rec;
    rec.push_back({0,0,0,0,0,0,0,0,0});

    long step = 0;
    while (t < t_max - rk_dt * 0.5) {
        rk4_step_linear(t, y, rk_dt, p);
        ++step;
        if (step % sample_interval == 0) {
            double ts = step / (double)sample_interval * sample_dt;
            Sample s;
            s.t = ts;
            s.zf = y[0]; s.vf = y[4];
            s.tf = y[1]; s.wf = y[5];
            s.zo = z_o_out(y, p);
            s.vo = zd_o_out(y, y + 4, p);
            s.to = y[3]; s.wo = y[7];
            rec.push_back(s);
        }
    }
    return rec;
}

// ====================================================================
// 【模块 L】XLSX 写入器 (STORE 方法, 无压缩)
// 功能: 手写 OOXML + ZIP (30 字节本地文件头 + 中央目录 + EOCD)
//       生成 Excel 可打开的 result3.xlsx
// ====================================================================
// 小端序写入 16 位/32 位无符号整数
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

// CRC32 校验和
uint32_t crc32_calc(const std::string& data) {
    uint32_t crc = 0xFFFFFFFF;
    for (unsigned char c : data) {
        crc ^= c;
        for (int k = 0; k < 8; ++k)
            crc = (crc >> 1) ^ (0xEDB88320 & -(crc & 1));
    }
    return ~crc;
}

struct ZipEntry {
    std::string name, data;
};

// 组装 ZIP 文件 (作为 xlsx 载体)
std::string make_xlsx(const std::vector<ZipEntry>& entries) {
    std::string out, cdir;
    std::vector<uint32_t> offsets;
    for (auto& e : entries) {
        uint32_t crc = crc32_calc(e.data);
        uint32_t usize = (uint32_t)e.data.size();
        offsets.push_back((uint32_t)out.size());

        // 本地文件头 (30 字节固定头)
        out += "PK\x03\x04";                    // 签名
        put_u16(out, 20);                       // 提取所需版本
        put_u16(out, 0);                        // 通用标志
        put_u16(out, 0);                        // 压缩方法 = STORE(0)
        put_u16(out, 0);                        // 修改时间
        put_u16(out, 0);                        // 修改日期
        put_u32(out, crc);                      // CRC32
        put_u32(out, usize);                    // 压缩后大小 (=未压缩)
        put_u32(out, usize);                    // 未压缩大小
        put_u16(out, (uint16_t)e.name.size());  // 文件名长度
        put_u16(out, 0);                        // 扩展字段长度
        out += e.name;
        out += e.data;                          // 数据 (未压缩)

        // 中央目录条目 (46 字节固定头)
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
        put_u32(cdir, offsets.back());          // 本地文件头偏移
        cdir += e.name;
    }

    uint32_t cdir_offset = (uint32_t)out.size();
    out += cdir;

    // EOCD (22 字节)
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

// 将仿真结果写入 result3.xlsx (9 列, 与结果需求 A~I 列序一致)
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
    // 第 1 行: 列标题 (inlineStr 类型)
    sheet << "<row r=\"1\">";
    for (int j = 0; j < 9; ++j)
        sheet << "<c r=\"" << cols[j] << "1\" t=\"inlineStr\"><is><t>"
              << headers[j] << "</t></is></c>";
    sheet << "</row>";

    // 第 2 行起: 数据 (纯数值)
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

// ====================================================================
// 【模块 M】CSV 结果写入 (UTF-8 BOM, 保证 Excel/VSCode 中文不乱码)
// ====================================================================
void write_result_csv(const std::string& path, const std::vector<Sample>& rec) {
    std::ofstream f(path, std::ios::binary);
    f << "\xEF\xBB\xBF";   // UTF-8 BOM, 使 Excel 正确识别中文
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

// 写 UTF-8 BOM
void put_bom(std::ofstream& out) { out << "\xEF\xBB\xBF"; }

// ====================================================================
// 【模块 N】路径解析
// ====================================================================
// 搜索附件文件所在目录 (支持多层目录结构)
std::string resolve_data_dir() {
    std::string c[] = {"题目信息/A题", "../../../题目信息/A题", "../题目信息/A题"};
    for (auto& s : c)
        if (std::filesystem::exists(s + "/附件3.xlsx")
            && std::filesystem::exists(s + "/附件4.xlsx")) return s;
    return "题目信息/A题";
}
// 自动递增运行输出目录 run001, run002, ...
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

// 在采样序列中找最接近目标时间 t_target 的索引
int find_idx(double t_target, const std::vector<Sample>& rec) {
    int best = 0;
    double best_diff = 1e30;
    for (int i = 0; i < (int)rec.size(); ++i) {
        double d = std::fabs(rec[i].t - t_target);
        if (d < best_diff) { best = i; best_diff = d; }
    }
    return best;
}

}  // namespace (匿名命名空间结束)

// ====================================================================
// main — 程序入口
// 命令行: --data <附件目录> --out <输出目录>
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

    // 控制台输出参数摘要
    std::cout << "======== Q3 参数 ========" << std::endl;
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

    // ==================== 正式积分 (L2 耦合模型) ====================
    std::cout << "\n======== 正式积分 (L2 耦合模型, dt=0.001 s) ========" << std::endl;
    FullRec fr = simulate_full(p, 0.001);
    std::cout << "输出点数: " << fr.samples.size() << " (含 t=0), 最后时刻 t="
              << fr.samples.back().t << " s" << std::endl;

    // ==================== result3.xlsx / result3.csv ====================
    write_result_xlsx(out_dir + "/result3.xlsx", fr.samples);
    std::cout << "\n已写入 " << out_dir << "/result3.xlsx" << std::endl;
    write_result_csv(out_dir + "/result3.csv", fr.samples);
    std::cout << "已写入 " << out_dir << "/result3.csv" << std::endl;

    // ==================== 论文五时刻 CSV ====================
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

        std::cout << "\n====== Q3 论文五时刻结果 ======" << std::endl;
        std::cout << "     t      zf     żf     θf    θ̇f    zo_out  żo_out    θo    θ̇o" << std::endl;
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

    // ==================== 弹簧长度检查 (硬约束 ξ>-0.2019575 m) ====================
    {
        double min_spring = 1e30, max_spring = -1e30;
        bool violated = false;
        for (size_t i = 0; i < fr.spring_length.size(); ++i) {
            if (fr.spring_length[i] < min_spring) min_spring = fr.spring_length[i];
            if (fr.spring_length[i] > max_spring) max_spring = fr.spring_length[i];
            if (fr.spring_length[i] <= p.H_o * 0.5) violated = true;
        }
        std::cout << "\n====== 弹簧长度检查（硬约束 ξ>-0.2019575 m）======" << std::endl;
        std::cout << "最小 ℓ=" << std::setprecision(10) << min_spring
                  << " m (下界=" << (p.H_o * 0.5) << " m)" << std::endl;
        std::cout << "ξ 最小值=" << (min_spring - p.l_e) << " m" << std::endl;
        if (violated) std::cout << "** 警告：弹簧几何长度越界！结果无效！ **" << std::endl;
        else std::cout << "弹簧长度全程有效。" << std::endl;

        std::ofstream f(out_dir + "/result3_弹簧长度检查.csv", std::ios::binary);
        put_bom(f);
        f << "时间t(s),弹簧瞬时长度l(m),滑移量xi(m)\n";
        for (size_t i = 0; i < fr.samples.size(); ++i)
            f << std::setprecision(10) << fr.samples[i].t << "," << fr.spring_length[i] << ","
              << (fr.spring_length[i] - p.l_e) << "\n";
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_弹簧长度检查.csv" << std::endl;
    }

    // ==================== 质量矩阵正定性检查 (硬约束 M 正定) ====================
    {
        double min_ev_all = 1e30;
        int neg_count = 0;
        for (size_t i = 0; i < fr.min_eigenvalue.size(); ++i) {
            if (fr.min_eigenvalue[i] < min_ev_all) min_ev_all = fr.min_eigenvalue[i];
            if (fr.min_eigenvalue[i] <= 0) ++neg_count;
        }
        std::cout << "\n====== 质量矩阵正定性检查（硬约束 min_eig(M)>0）======" << std::endl;
        std::cout << "全程最小特征值=" << std::setprecision(10) << min_ev_all << std::endl;
        if (neg_count > 0) std::cout << "** 警告：" << neg_count << " 个输出时刻 M 非正定！ **" << std::endl;
        else std::cout << "所有输出时刻 M 正定。" << std::endl;

        std::ofstream f(out_dir + "/result3_质量矩阵正定性.csv", std::ios::binary);
        put_bom(f);
        f << "时间t(s),质量矩阵最小特征值min_eig_M\n";
        for (size_t i = 0; i < fr.samples.size(); ++i)
            f << std::setprecision(10) << fr.samples[i].t << "," << fr.min_eigenvalue[i] << "\n";
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_质量矩阵正定性.csv" << std::endl;
    }

    // ==================== 能量恒等式检查 (式 Q3-18) ====================
    {
        double max_relerr = 0.0;
        for (size_t i = 0; i < fr.energy.size(); ++i) {
            double E = fr.energy[i];
            double bal = fr.energy_balance[i];
            double relerr = std::fabs(E - bal) / (1.0 + std::fabs(bal));
            if (relerr > max_relerr) max_relerr = relerr;
        }
        std::cout << "\n====== 能量恒等式检查（式 Q3-18）======" << std::endl;
        std::cout << "max |E(t) - [E(0)+∫Ė dt]| / (1+|E(0)+∫Ė dt|) = "
                  << std::setprecision(10) << max_relerr << std::endl;
        if (max_relerr <= 1e-3)
            std::cout << "通过 (≤ 1e-3)。" << std::endl;
        else
            std::cout << "未通过 (> 1e-3)。" << std::endl;

        std::ofstream f(out_dir + "/result3_能量恒等式.csv", std::ios::binary);
        put_bom(f);
        f << "时间t(s),机械能E(J),能量对账E0+∫Ėdt(J),相对残差\n";
        for (size_t i = 0; i < fr.energy.size(); ++i) {
            double relerr = std::fabs(fr.energy[i] - fr.energy_balance[i])
                          / (1.0 + std::fabs(fr.energy_balance[i]));
            f << std::setprecision(10) << fr.samples[i].t << "," << fr.energy[i] << ","
              << fr.energy_balance[i] << "," << relerr << "\n";
        }
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_能量恒等式.csv" << std::endl;
    }

    // ==================== 线性模型对照 (诊断, 不覆盖 result3.xlsx) ====================
    std::cout << "\n======== 线性模型对照 (诊断, dt=0.001 s) ========" << std::endl;
    auto lin_rec = simulate_linear(p, 0.001);
    {
        std::ofstream f(out_dir + "/result3_线性模型对照.csv", std::ios::binary);
        put_bom(f);
        f << "时刻t(s),浮子垂荡位移差dzf(m),浮子垂荡速度差(m/s),浮子纵摇角差dθf(rad),"
          << "浮子纵摇角速度差(rad/s),振子垂荡位移差dzo(m),振子垂荡速度差(m/s),"
          << "振子纵摇角差dθo(rad),振子纵摇角速度差(rad/s)\n";
        for (double kt : key_times) {
            int ni = find_idx(kt, fr.samples);
            int li = find_idx(kt, lin_rec);
            f << kt << ","
              << (fr.samples[ni].zf - lin_rec[li].zf) << ","
              << (fr.samples[ni].vf - lin_rec[li].vf) << ","
              << (fr.samples[ni].tf - lin_rec[li].tf) << ","
              << (fr.samples[ni].wf - lin_rec[li].wf) << ","
              << (fr.samples[ni].zo - lin_rec[li].zo) << ","
              << (fr.samples[ni].vo - lin_rec[li].vo) << ","
              << (fr.samples[ni].to - lin_rec[li].to) << ","
              << (fr.samples[ni].wo - lin_rec[li].wo) << "\n";
        }
        f.close();
        std::cout << "已写入 " << out_dir << "/result3_线性模型对照.csv" << std::endl;

        double max_zf_diff = 0, max_tf_diff = 0, max_zo_diff = 0, max_to_diff = 0;
        for (size_t i = 0; i < fr.samples.size() && i < lin_rec.size(); ++i) {
            double d1 = std::fabs(fr.samples[i].zf - lin_rec[i].zf);
            double d2 = std::fabs(fr.samples[i].tf - lin_rec[i].tf);
            double d3 = std::fabs(fr.samples[i].zo - lin_rec[i].zo);
            double d4 = std::fabs(fr.samples[i].to - lin_rec[i].to);
            if (d1 > max_zf_diff) max_zf_diff = d1;
            if (d2 > max_tf_diff) max_tf_diff = d2;
            if (d3 > max_zo_diff) max_zo_diff = d3;
            if (d4 > max_to_diff) max_to_diff = d4;
        }
        std::cout << "最大差值 |Δzf|=" << std::setprecision(6) << max_zf_diff
                  << " m  |Δθf|=" << max_tf_diff << " rad  |Δzo|=" << max_zo_diff
                  << " m  |Δθo|=" << max_to_diff << " rad" << std::endl;
    }

    // ==================== 响应幅值摘要 ====================
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

    std::cout << "\n======== Q3 求解完成 ========" << std::endl;
    std::cout << "全部结果输出至: " << out_dir << std::endl;

    return 0;
}
