// Q1 v1 求解器（C++ 详细注释版）. Owner: 罗懿 (L, 建模手).
//
// 数据来源（硬性要求：全部来自真实附件文件，杜绝硬编码捏造）：
//   题目信息/A题/附件3.xlsx —— 在频率 ω=1.4005 行读取：垂荡附加质量 m_a、垂荡兴波阻尼系数 b_h、垂子激励力振幅 F
//   题目信息/A题/附件4.xlsx —— 读取：浮子质量 m_f、振子质量 m_o、浮子底半径(水线面半径) R_f、海水密度 ρ、重力加速度 g、弹簧刚度 k
// 模型依据：模型/Q1/v1/数学模型.md 式(Q1-1)~(Q1-6)。
// 统一耗散函数 D(v_r)=c·v_r（情形一）或 D(v_r)=c0·|v_r|^p·v_r（情形二），两种情形共用同一 RK4 积分器。
// 题面场景常量（非文件测量数据，来自题面/约束.md）：ω=1.4005 s⁻¹，情形一 c=10000 N·s/m，情形二 c0=10000、p=0.5。

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
const double OMEGA = 1.4005;   // 题面 Q1 场景：入射波圆频率（s⁻¹）
const double C0 = 10000.0;     // 题面 Q1 情形一常量阻尼系数 c / 情形二幂律比例系数 c0（N·s/m 或相应量纲）
const double P = 0.5;          // 题面 Q1 情形二幂指数 p

// 物理参数（运行期由附件文件读取填入，代码内不写死数值）
struct Params {
    double m_f, m_o, m_a, b_h, F, rho, g, R_f, k;
};

// 一行采样结果：时间、浮子位移、浮子速度、振子位移、振子速度
struct Sample {
    double t, zf, vf, zo, vo;
};

// 读取整个文件为二进制字符串
std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("无法打开文件: " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// 对 ZIP 内 raw deflate 流解压（xlsx 内 sheet/sharedStrings 通常为 deflate 压缩）
std::string inflate_raw(const std::string& in) {
    z_stream zs;
    zs.zalloc = Z_NULL; zs.zfree = Z_NULL; zs.opaque = Z_NULL;
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

// 一个单元格：区分字符串（共享字符串索引或内联串）与数值
struct Cell {
    bool is_str = false;
    double num = 0.0;
    std::string str;
};

// 解析 sharedStrings.xml 为共享字符串表（合并同一 <si> 内的多段 <t> 文本）
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

// 解析 sheet1.xml：提取每个单元格引用(如 "B2") 到 Cell 的映射
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
        cell.is_str = (c.find(" t=\"s\"") != std::string::npos) || (c.find(" t=\"inlineStr\"") != std::string::npos);
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

// 【重要模块】xlsx(ZIP)解析：枚举中央目录条目，按需解压 sheet1.xml 与 sharedStrings.xml。
//   严格只读取题面给定的真实附件，不读取无关文件，且文件名精确匹配。
std::map<std::string, std::string> unzip_entries(const std::string& data) {
    std::map<std::string, std::string> out;
    size_t eocd = std::string::npos;
    for (size_t i = data.size() >= 22 ? data.size() - 22 : 0; ; --i) {
        if (i + 4 <= data.size() &&
            (uint8_t)data[i] == 0x50 && (uint8_t)data[i + 1] == 0x4b &&
            (uint8_t)data[i + 2] == 0x05 && (uint8_t)data[i + 3] == 0x06) {
            eocd = i;
            break;
        }
        if (i == 0) break;
    }
    if (eocd == std::string::npos) throw std::runtime_error("不是有效的 xlsx(zip) 文件");
    uint32_t central_off = (uint32_t)(uint8_t)data[eocd + 16]
        | ((uint32_t)(uint8_t)data[eocd + 17] << 8)
        | ((uint32_t)(uint8_t)data[eocd + 18] << 16)
        | ((uint32_t)(uint8_t)data[eocd + 19] << 24);
    size_t p = central_off;
    while (p + 4 <= data.size() &&
           (uint8_t)data[p] == 0x50 && (uint8_t)data[p + 1] == 0x4b &&
           (uint8_t)data[p + 2] == 0x01 && (uint8_t)data[p + 3] == 0x02) {
        uint16_t method = (uint16_t)((uint8_t)data[p + 10] | ((uint8_t)data[p + 11] << 8));
        uint32_t comp_size = (uint32_t)(uint8_t)data[p + 20]
            | ((uint32_t)(uint8_t)data[p + 21] << 8)
            | ((uint32_t)(uint8_t)data[p + 22] << 16)
            | ((uint32_t)(uint8_t)data[p + 23] << 24);
        uint16_t fnlen = (uint16_t)((uint8_t)data[p + 28] | ((uint8_t)data[p + 29] << 8));
        uint16_t extralen = (uint16_t)((uint8_t)data[p + 30] | ((uint8_t)data[p + 31] << 8));
        uint16_t commentlen = (uint16_t)((uint8_t)data[p + 32] | ((uint8_t)data[p + 33] << 8));
        uint32_t local_off = (uint32_t)(uint8_t)data[p + 42]
            | ((uint32_t)(uint8_t)data[p + 43] << 8)
            | ((uint32_t)(uint8_t)data[p + 44] << 16)
            | ((uint32_t)(uint8_t)data[p + 45] << 24);
        std::string name = data.substr(p + 46, fnlen);
        size_t lp = local_off;
        uint16_t lfnlen = (uint16_t)((uint8_t)data[lp + 26] | ((uint8_t)data[lp + 27] << 8));
        uint16_t lextralen = (uint16_t)((uint8_t)data[lp + 28] | ((uint8_t)data[lp + 29] << 8));
        size_t datastart = lp + 30 + lfnlen + lextralen;
        std::string comp = data.substr(datastart, comp_size);
        out[name] = (method == 0) ? comp : inflate_raw(comp);
        p += 46 + fnlen + extralen + commentlen;
    }
    return out;
}

// 读取 xlsx 返回单元格映射；优先使用 sheet1.xml
std::map<std::string, Cell> read_xlsx(const std::string& path) {
    std::string data = read_file(path);
    auto entries = unzip_entries(data);
    auto sit = entries.find("xl/sharedStrings.xml");
    std::vector<std::string> ss;
    if (sit != entries.end()) ss = parse_shared_strings(sit->second);
    auto it = entries.find("xl/worksheets/sheet1.xml");
    if (it == entries.end()) it = entries.find("xl/worksheets/sheet1.xml");
    if (it == entries.end()) throw std::runtime_error("xlsx 缺少 sheet1.xml: " + path);
    return parse_sheet(it->second, ss);
}

// 取单元格列的字母部分（如 "B2" -> "B"）
std::string col_letters(const std::string& ref) {
    size_t i = 0;
    while (i < ref.size() && std::isalpha((unsigned char)ref[i])) ++i;
    return ref.substr(0, i);
}
// 取单元格行号（如 "B2" -> 2）
int row_number(const std::string& ref) {
    size_t i = 0;
    while (i < ref.size() && std::isalpha((unsigned char)ref[i])) ++i;
    return std::stoi(ref.substr(i));
}
// 归一化参数标签：去掉末尾 " (单位)" 后缀，使其与模型符号名一致（附件标签含单位，如 "浮子质量 (kg)"）
std::string norm_label(const std::string& s) {
    size_t pos = s.rfind(" (");
    if (pos != std::string::npos) return s.substr(0, pos);
    return s;
}

// 【重要模块】参数读取（数据来源审计）：从附件4读取质量/几何/密度/重力/弹簧，从附件3按频率 ω 定位 Q1 行读取水动力参数。
//   对每个读取项做存在性与有限值断言；附件3 频率行必须与题面 ω 精确匹配，否则报错（防止读错文件/错行）。
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
        if (!std::isfinite(it->second)) throw std::runtime_error("附件4 参数非有限: " + name);
        return it->second;
    };
    Params p;
    p.m_f = need4("浮子质量");
    p.R_f = need4("浮子底半径");
    p.m_o = need4("振子质量");
    p.rho = need4("海水的密度");
    p.g = need4("重力加速度");
    p.k = need4("弹簧刚度");

    auto c3 = read_xlsx(attach3);
    bool found = false;
    for (auto& kv : c3) {
        if (col_letters(kv.first) == "B" && !kv.second.is_str) {
            int r = row_number(kv.first);
            if (std::fabs(kv.second.num - OMEGA) < 1e-6) {
                auto get3 = [&](const std::string& col) -> double {
                    std::string ref = col + std::to_string(r);
                    auto it = c3.find(ref);
                    if (it == c3.end() || it->second.is_str)
                        throw std::runtime_error("附件3 缺少列 " + col + " 行 " + std::to_string(r));
                    if (!std::isfinite(it->second.num))
                        throw std::runtime_error("附件3 参数非有限 列 " + col);
                    return it->second.num;
                };
                p.m_a = get3("C");   // 垂荡附加质量 m_a
                p.b_h = get3("E");   // 垂荡兴波阻尼系数 b_h
                p.F = get3("G");     // 垂荡激励力振幅 F
                found = true;
                break;
            }
        }
    }
    if (!found) throw std::runtime_error("附件3 未找到频率 1.4005 对应的 Q1 行");
    return p;
}

// 【重要模块】阻尼/PTO 耗散函数：负责把相对速度映射为 PTO 阻尼力标量 D(v_r)。
//   情形一 D=c·v_r（线性）；情形二 D=c0·|v_r|^p·v_r（幂律）。|v_r|^p·v_r 与 v_r 同号 => D·v_r≥0（耗散非负，绝不向系统注入能量）。
double damping(double vr, int type) {
    if (type == 1) return C0 * vr;
    return C0 * std::pow(std::fabs(vr), P) * vr;
}

// 【重要模块】运动方程右端项：负责把 数学模型.md 式(Q1-2)(Q1-3) 写成一阶状态导数 dx/dt。
//   状态 x=(z_f, ż_f, z_o, ż_o)；相对位移 ξ=z_o−z_f，相对速度 v_r=ż_o−ż_f；静水恢复刚度 K_h=ρ g π R_f²（已在式中内联）。
void rhs(double t, const double x[4], double dx[4], int type, const Params& p) {
    double xi = x[2] - x[0];          // 振子相对浮子位移 ξ
    double vr = x[3] - x[1];          // 相对速度 v_r
    double D = damping(vr, type);     // PTO 阻尼力标量
    double exc = p.F * std::cos(OMEGA * t);  // 波浪激励力 f cos(ωt)
    dx[0] = x[1];
    dx[1] = (exc - p.b_h * x[1] - (p.rho * p.g * PI * p.R_f * p.R_f) * x[0] + p.k * xi + D) / (p.m_f + p.m_a);
    dx[2] = x[3];
    dx[3] = (-p.k * xi - D) / p.m_o;
}

// 【重要模块】经典四阶 Runge-Kutta 积分器：对一阶状态方程单步推进，全局误差 O(dt⁴)。
void rk4(double t, double x[4], double dt, int type, const Params& p) {
    double k1[4], k2[4], k3[4], k4[4], xt[4];
    rhs(t, x, k1, type, p);
    for (int i = 0; i < 4; ++i) xt[i] = x[i] + 0.5 * dt * k1[i];
    rhs(t + 0.5 * dt, xt, k2, type, p);
    for (int i = 0; i < 4; ++i) xt[i] = x[i] + 0.5 * dt * k2[i];
    rhs(t + 0.5 * dt, xt, k3, type, p);
    for (int i = 0; i < 4; ++i) xt[i] = x[i] + dt * k3[i];
    rhs(t + dt, xt, k4, type, p);
    for (int i = 0; i < 4; ++i) x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

// 【重要模块】时间推进与采样：积分前 40 个波浪周期，按题面 0.2 s 间隔采样（dt=0.001 远细于采样间隔保证精度）。
std::vector<Sample> simulate(int type, const Params& p, double dt) {
    int steps_per_sample = (int)std::lround(0.2 / dt);
    const double T = 2.0 * PI / OMEGA;
    const double tmax = 40.0 * T;
    double x[4] = {0.0, 0.0, 0.0, 0.0};  // 静水平衡初值：位移、速度全为零
    std::vector<Sample> out;
    long n = 0;
    double t = 0.0;
    out.push_back({0.0, x[0], x[1], x[2], x[3]});  // 记录 t=0
    while (t < tmax - 1e-12) {
        rk4(t, x, dt, type, p);
        t += dt;
        ++n;
        if (n % steps_per_sample == 0) {            // 到 0.2 s 倍数时采样
            double ts = (n / steps_per_sample) * 0.2;
            out.push_back({ts, x[0], x[1], x[2], x[3]});
        }
    }
    return out;
}

uint32_t crc32(const std::string& data) {
    uint32_t crc = 0xFFFFFFFFu;
    for (unsigned char c : data) {
        crc ^= c;
        for (int i = 0; i < 8; ++i) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
    }
    return ~crc;
}
void put16(std::ostream& os, uint16_t v) {
    os.put((char)(v & 0xFF));
    os.put((char)((v >> 8) & 0xFF));
}
void put32(std::ostream& os, uint32_t v) {
    os.put((char)(v & 0xFF));
    os.put((char)((v >> 8) & 0xFF));
    os.put((char)((v >> 16) & 0xFF));
    os.put((char)((v >> 24) & 0xFF));
}
struct ZipEntry { std::string name; std::string data; };

// 【重要模块】xlsx 写出：把采样结果写成真实 Office Open XML（.xlsx）工作簿，表头与题面模板一致并保留合并单元格。
void write_xlsx(const std::string& filename, const std::vector<Sample>& s) {
    const std::string head =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
        "<sheetData>"
        "<row r=\"1\">"
        "<c r=\"A1\" t=\"inlineStr\"><is><t>时间 (s)</t></is></c>"
        "<c r=\"B1\" t=\"inlineStr\"><is><t>浮子</t></is></c>"
        "<c r=\"C1\"/>"
        "<c r=\"D1\" t=\"inlineStr\"><is><t>振子</t></is></c>"
        "<c r=\"E1\"/>"
        "</row>"
        "<row r=\"2\">"
        "<c r=\"A2\"/>"
        "<c r=\"B2\" t=\"inlineStr\"><is><t>位移 (m)</t></is></c>"
        "<c r=\"C2\" t=\"inlineStr\"><is><t>速度 (m/s)</t></is></c>"
        "<c r=\"D2\" t=\"inlineStr\"><is><t>位移 (m)</t></is></c>"
        "<c r=\"E2\" t=\"inlineStr\"><is><t>速度 (m/s)</t></is></c>"
        "</row>";
    std::ostringstream body;
    int r = 3;
    for (const auto& p : s) {
        std::ostringstream v;
        v << std::setprecision(12) << p.t << "</v></c>"
          << "<c r=\"B" << r << "\"><v>" << std::setprecision(12) << p.zf << "</v></c>"
          << "<c r=\"C" << r << "\"><v>" << std::setprecision(12) << p.vf << "</v></c>"
          << "<c r=\"D" << r << "\"><v>" << std::setprecision(12) << p.zo << "</v></c>"
          << "<c r=\"E" << r << "\"><v>" << std::setprecision(12) << p.vo << "</v></c>";
        body << "<row r=\"" << r << "\">"
             << "<c r=\"A" << r << "\"><v>" << v.str() << "</row>";
        ++r;
    }
    const std::string sheet_full = head + body.str() + "</sheetData>"
        "<mergeCells count=\"3\">"
        "<mergeCell ref=\"B1:C1\"/><mergeCell ref=\"D1:E1\"/><mergeCell ref=\"A1:A2\"/>"
        "</mergeCells></worksheet>";

    std::vector<ZipEntry> entries = {
        {"[Content_Types].xml",
         "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
         "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
         "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
         "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
         "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
         "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
         "<Override PartName=\"/xl/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>"
         "</Types>"},
        {"_rels/.rels",
         "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
         "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
         "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
         "</Relationships>"},
        {"xl/workbook.xml",
         "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
         "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
         "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">"
         "<sheets><sheet name=\"Sheet1\" sheetId=\"1\" r:id=\"rId1\"/></sheets></workbook>"},
        {"xl/_rels/workbook.xml.rels",
         "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
         "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
         "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
         "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>"
         "</Relationships>"},
        {"xl/styles.xml",
         "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
         "<styleSheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
         "<fonts count=\"1\"><font><sz val=\"11\"/><name val=\"Calibri\"/></font></fonts>"
         "<fills count=\"1\"><fill><patternFill patternType=\"none\"/></fill></fills>"
         "<borders count=\"1\"><border/></borders>"
         "<cellStyleXfs count=\"1\"><xf/></cellStyleXfs>"
         "<cellXfs count=\"1\"><xf/></cellXfs>"
         "</styleSheet>"},
        {"xl/worksheets/sheet1.xml", sheet_full},
    };
    std::ofstream out(filename, std::ios::binary);
    std::vector<uint32_t> offsets, crcs;
    for (const auto& e : entries) {
        uint32_t crc = crc32(e.data);
        uint32_t offset = (uint32_t)out.tellp();
        out.put(0x50); out.put(0x4b); out.put(0x03); out.put(0x04);
        put16(out, 20); put16(out, 0); put16(out, 0); put16(out, 0); put16(out, 0);
        put32(out, crc);
        put32(out, (uint32_t)e.data.size());
        put32(out, (uint32_t)e.data.size());
        put16(out, (uint16_t)e.name.size());
        put16(out, 0);
        out.write(e.name.data(), e.name.size());
        out.write(e.data.data(), e.data.size());
        offsets.push_back(offset); crcs.push_back(crc);
    }
    uint32_t central_offset = (uint32_t)out.tellp();
    for (size_t i = 0; i < entries.size(); ++i) {
        const auto& e = entries[i];
        out.put(0x50); out.put(0x4b); out.put(0x01); out.put(0x02);
        put16(out, 20); put16(out, 20); put16(out, 0); put16(out, 0); put16(out, 0); put16(out, 0);
        put32(out, crcs[i]);
        put32(out, (uint32_t)e.data.size());
        put32(out, (uint32_t)e.data.size());
        put16(out, (uint16_t)e.name.size());
        put16(out, 0); put16(out, 0); put16(out, 0); put16(out, 0); put32(out, 0);
        put32(out, offsets[i]);
        out.write(e.name.data(), e.name.size());
    }
    uint32_t central_size = (uint32_t)out.tellp() - central_offset;
    out.put(0x50); out.put(0x4b); out.put(0x05); out.put(0x06);
    put16(out, 0); put16(out, 0);
    put16(out, (uint16_t)entries.size());
    put16(out, (uint16_t)entries.size());
    put32(out, central_size); put32(out, central_offset); put16(out, 0);
    out.close();
}

// 【重要模块】关键时间点 CSV（utf-8-sig，含中文表头与单位）：供论文直接引用 10/20/40/60/100 s 的位移与速度。
void write_key_csv(const std::string& filename, const std::vector<Sample>& s) {
    double times[5] = {10.0, 20.0, 40.0, 60.0, 100.0};
    std::ofstream out(filename, std::ios::binary);
    out << "\xEF\xBB\xBF";
    out << "时间(s),浮子位移(m),浮子速度(m/s),振子位移(m),振子速度(m/s)\n";
    for (double tt : times) {
        long idx = (long)std::llround(tt / 0.2);
        const Sample& q = s[idx];
        out << std::setprecision(12) << q.t << "," << q.zf << "," << q.vf << "," << q.zo << "," << q.vo << "\n";
    }
    out.close();
}

// 【重要模块】PTO 能量输出系统功率与累计能量 CSV：瞬时功率 P=D·v_r，累计能量用梯形法数值积分。
//   作为 Q1 物理补充结果，也为 Q2（最大平均输出功率）提供前置量。
void write_power_csv(const std::string& filename, int type, const std::vector<Sample>& s) {
    std::ofstream out(filename, std::ios::binary);
    out << "\xEF\xBB\xBF";
    out << "时间(s),PTO能量输出系统瞬时功率(瓦特),PTO能量输出系统累计能量(焦耳)\n";
    double E = 0.0;
    double prevP = 0.0;
    for (size_t i = 0; i < s.size(); ++i) {
        double vr = s[i].vo - s[i].vf;
        double P = damping(vr, type) * vr;
        if (i > 0) E += 0.5 * (P + prevP) * (s[i].t - s[i - 1].t);
        out << std::setprecision(12) << s[i].t << "," << P << "," << E << "\n";
        prevP = P;
    }
    out.close();
}

// 【重要模块】两情形对比 CSV：同采样时刻 case1 减 case2 的位移/速度差，用于灵敏度/对比分析。
void write_compare_csv(const std::string& filename, const std::vector<Sample>& s1, const std::vector<Sample>& s2) {
    std::ofstream out(filename, std::ios::binary);
    out << "\xEF\xBB\xBF";
    out << "时间(s),浮子位移差_情形1减情形2(m),浮子速度差_情形1减情形2(m/s),振子位移差_情形1减情形2(m),振子速度差_情形1减情形2(m/s)\n";
    size_t n = s1.size() < s2.size() ? s1.size() : s2.size();
    for (size_t i = 0; i < n; ++i) {
        out << std::setprecision(12) << s1[i].t << ","
            << (s1[i].zf - s2[i].zf) << "," << (s1[i].vf - s2[i].vf) << ","
            << (s1[i].zo - s2[i].zo) << "," << (s1[i].vo - s2[i].vo) << "\n";
    }
    out.close();
}

// 【重要模块】数值诊断 CSV（验证/误差分析）：幅度统计、耗散功率非负断言、RK4 步长收敛误差、静水恢复刚度。
void write_diag_csv(const std::string& filename, const Params& p, const std::vector<Sample>& s1,
                    const std::vector<Sample>& s2, const std::vector<Sample>& s1f, const std::vector<Sample>& s2f) {
    std::ofstream out(filename, std::ios::binary);
    out << "\xEF\xBB\xBF";
    out << "检查项,情形,数值\n";
    auto maxabs = [](const std::vector<Sample>& s, double Sample::*f) {
        double m = 0.0;
        for (const auto& q : s) m = std::max(m, std::fabs(q.*f));
        return m;
    };
    auto mindiss = [&](const std::vector<Sample>& s, int type) {
        double m = 1e300;
        for (const auto& q : s) {
            double vr = q.vo - q.vf;
            m = std::min(m, damping(vr, type) * vr);
        }
        return m;
    };
    out << "最大|浮子位移|(m),情形一," << std::setprecision(12) << maxabs(s1, &Sample::zf) << "\n";
    out << "最大|浮子速度|(m/s),情形一," << maxabs(s1, &Sample::vf) << "\n";
    out << "最大|振子位移|(m),情形一," << maxabs(s1, &Sample::zo) << "\n";
    out << "最大|振子速度|(m/s),情形一," << maxabs(s1, &Sample::vo) << "\n";
    out << "最大|浮子位移|(m),情形二," << maxabs(s2, &Sample::zf) << "\n";
    out << "最大|浮子速度|(m/s),情形二," << maxabs(s2, &Sample::vf) << "\n";
    out << "最大|振子位移|(m),情形二," << maxabs(s2, &Sample::zo) << "\n";
    out << "最大|振子速度|(m/s),情形二," << maxabs(s2, &Sample::vo) << "\n";
    out << "最小耗散功率D*vr(瓦特,应>=0),情形一," << mindiss(s1, 1) << "\n";
    out << "最小耗散功率D*vr(瓦特,应>=0),情形二," << mindiss(s2, 2) << "\n";
    double conv = 0.0;
    size_t n = s1.size() < s1f.size() ? s1.size() : s1f.size();
    for (size_t i = 0; i < n; ++i)
        conv = std::max(conv, std::fabs(s1[i].zf - s1f[i].zf));
    out << "RK4步长收敛误差_最大浮子位移差_粗步长0.001相对细步长0.0005(m),情形一," << conv << "\n";
    double conv2 = 0.0;
    n = s2.size() < s2f.size() ? s2.size() : s2f.size();
    for (size_t i = 0; i < n; ++i)
        conv2 = std::max(conv2, std::fabs(s2[i].zf - s2f[i].zf));
    out << "RK4步长收敛误差_最大浮子位移差_粗步长0.001相对细步长0.0005(m),情形二," << conv2 << "\n";
    out << "静水恢复刚度Kh(N/m),-, " << (p.rho * p.g * PI * p.R_f * p.R_f) << "\n";
    out.close();
    if (mindiss(s1, 1) < -1e-9) throw std::runtime_error("情形一耗散功率为负，阻尼符号错误");
    if (mindiss(s2, 2) < -1e-9) throw std::runtime_error("情形二耗散功率为负，阻尼符号错误");
}

// 在控制台打印论文所需的关键时间点响应，便于核对（正式记录以 CSV 文件为准）。
void report_paper(const std::vector<Sample>& s, const char* label) {
    std::cout << "===== " << label << " =====" << std::endl;
    double times[5] = {10.0, 20.0, 40.0, 60.0, 100.0};
    for (double tt : times) {
        long idx = (long)std::llround(tt / 0.2);
        const Sample& q = s[idx];
        std::cout << std::setprecision(12)
                  << "t=" << q.t << " s  浮子位移=" << q.zf << " m  浮子速度=" << q.vf
                  << " m/s  振子位移=" << q.zo << " m  振子速度=" << q.vo << " m/s\n";
    }
}

// 解析附件目录候选（兼容从不同工作目录运行），返回第一个存在的路径
std::string resolve_data_dir() {
    std::string candidates[] = {"题目信息/A题", "../../../题目信息/A题", "../题目信息/A题"};
    for (const auto& c : candidates) {
        if (std::filesystem::exists(c + "/附件3.xlsx") && std::filesystem::exists(c + "/附件4.xlsx"))
            return c;
    }
    return "题目信息/A题";
}
// 解析输出目录：默认写入本版本对应的结果归档目录，自动避开已存在的 run 编号
std::string resolve_out_dir() {
    const std::string base = "结果（包括各种过程数据与审查）/运行输出";
    if (!std::filesystem::exists(base)) return base + "/run001";
    for (int n = 1; n <= 999; ++n) {
        std::string d = base + "/run" + (n < 10 ? "00" : (n < 100 ? "0" : "")) + std::to_string(n);
        if (!std::filesystem::exists(d)) return d;
    }
    return base + "/run999";
}

}  // namespace

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

    std::string attach3 = data_dir + "/附件3.xlsx";
    std::string attach4 = data_dir + "/附件4.xlsx";
    Params p = load_params(attach3, attach4);

    std::cout << "数据来源核对（均读取自附件文件，非硬编码）：" << std::endl;
    std::cout << "  附件4: 浮子质量=" << p.m_f << " 振子质量=" << p.m_o << " 浮子半径=" << p.R_f
              << " 密度=" << p.rho << " 重力=" << p.g << " 弹簧刚度=" << p.k << std::endl;
    std::cout << "  附件3(ω=" << OMEGA << "): 附加质量=" << p.m_a << " 兴波阻尼=" << p.b_h
              << " 激励力振幅=" << p.F << std::endl;

    auto s1 = simulate(1, p, 0.001);
    auto s2 = simulate(2, p, 0.001);
    auto s1f = simulate(1, p, 0.0005);
    auto s2f = simulate(2, p, 0.0005);

    write_xlsx(out_dir + "/result1-1.xlsx", s1);
    write_xlsx(out_dir + "/result1-2.xlsx", s2);
    write_key_csv(out_dir + "/result1-1_关键时间点.csv", s1);
    write_key_csv(out_dir + "/result1-2_关键时间点.csv", s2);
    write_power_csv(out_dir + "/result1-1_功率与能量.csv", 1, s1);
    write_power_csv(out_dir + "/result1-2_功率与能量.csv", 2, s2);
    write_compare_csv(out_dir + "/result1_两情形对比.csv", s1, s2);
    write_diag_csv(out_dir + "/result1_数值诊断.csv", p, s1, s2, s1f, s2f);

    report_paper(s1, "情形一 常量阻尼 c=10000 N·s/m");
    report_paper(s2, "情形二 幂律阻尼 c0=10000,p=0.5");
    std::cout << "采样点数: 情形一 " << s1.size() << " 行, 情形二 " << s2.size() << " 行 (含 t=0)" << std::endl;
    std::cout << "输出目录: " << out_dir << std::endl;
    return 0;
}
