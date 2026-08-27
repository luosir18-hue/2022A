# Stage QA: paper-model-distillation

owner_name: 罗懿
owner_code: L
owner_role: 建模手
version: v2
status: blocked
based_on: 模型/Q2/v2/论文数学模型.md；模型/Q2/v2/数学模型.md；用户于 2026-08-26 给出的论文表达边界

AI status: blocked
Human spot-check: pending
Stage artifacts checked:
- 模型/Q2/v2/论文数学模型.md
- 模型/Q2/v2/模型简报.md#论文数学模型入口
- 模型/Q2/v2/数学模型.md
- 模型/Q2/v2/建模思路.md
- 模型/Q2/v2/结果需求.md
- 模型/Q2/v2/约束.md

## Scorecard

| Item | AI score 0/1/2 | Evidence | AI fix if not 2 | Modeler spot-check |
|---|---:|---|---|---|
| Same-version source | 2 | 派生稿登记 Q2/v2 五个模型源文件及 SHA-256 | — | pending |
| Forecast score coverage | 0 | `整体思路_overall_idea.md#写作得分点预测` 尚无可识别 `score_id` | problem-intake 补齐后逐项回填覆盖表 | pending |
| Problem/model distinction | 2 | 正文分列“题面给定条件”和“模型必要假设与定义” | — | pending |
| Model closure | 2 | 变量、目标、硬约束、初值、可行域、跨问依赖和输出均有正文锚点 | — | pending |
| Formula/constraint coverage | 2 | 内部覆盖表逐项映射源模型 Q2-1 至 Q2-11 | — | pending |
| Meaning preservation | 2 | 常量阻尼实数简谐表达与详细模型一元功率函数等价；非线性目标仍取 RK4 时域积分 | — | pending |
| Concise paper expression | 2 | 正文按 Q1 动力学—功率—平均—解析求导—RK4 搜索组织 | — | pending |
| Formula semantic lead-in | 2 | 机械检查未再报告公式前置语义问题 | — | pending |
| Parameterized equations | 2 | 控制方程保持参数形式，题面数值集中在题面条件段 | — | pending |
| Selective numbering | 2 | 核心方程、功率、目标、更新和输出编号；局部定义与代数中间式不编号 | — | pending |
| Formula readability | 2 | 未使用复阻抗、复振幅、Gamma 函数、谐波平衡或等效线性化；导数使用 $d/dt$ | — | pending |
| Clean paper prose | 2 | 正文未检出修改痕迹、程序保护旁白和内部审计语言 | — | pending |
| Omission traceability | 2 | 复振幅表达、单频辅助核对和代码级参数均在内部登记去向 | — | pending |
| Writer entry | 2 | `模型简报.md#论文数学模型入口` 已登记路径、读取顺序和 stale 条件 | — | pending |
| Human semantic spot-check | 1 | 用户已明确整体保留/省略方向，但尚未审阅本次落盘正文 | 用户核对实数简谐推导及正文边界 | pending |

AI total: 27/30
AI threshold: blocked because one required item scores 0
Human gate status: pending

## Mathematical spot-check

- 在 $c=0,10000,37193.8119,100000$ 四个点上，将正文实数简谐振幅公式与原二自由度四元实系数方程直接求解进行比较；相对速度振幅误差不超过 $2.3\times10^{-16}$，平均功率误差不超过 $5.5\times10^{-13}\ \mathrm{W}$。
- 由 $d\overline P_1/dc$ 的符号可知平均功率在 $c_s=\sqrt{A/C}$ 左侧递增、右侧递减；代入本问参数得到 $c_s=37193.8119$、$P_1^*=229.3339\ \mathrm{W}$。
- L-modeling 机械检查当前只剩“写作得分点预测无 `score_id`”这一阻塞项；公式语义、正文结构、条件分层与模型闭环的机械警告已清零。

## Red Flags

- 写作得分点预测尚未建立，无法完成官方论文位置的逐项覆盖审计。
- `结果需求.md#论文输出要求`、算法步长与收敛容差尚未填写，Q2/v2 也没有同版本正式运行证据，因此正文不报告速度相关阻尼的最优数值。
- Q1/Q2 源模型中的关键假设仍处于人工确认门禁；当前稿不能标记为正式采用版本。

## AI Suggestions

1. 建议先由模型员逐段核对“常量阻尼实数推导”和“非线性阻尼 RK4 搜索”是否与答辩口径一致。依据是本次改写改变了论文呈现方式但未改变数学关系；风险是未经人工抽查便进入论文会违反版本门禁；下一检查是将 `human_spot_check` 更新为 `pass` 或列出需改句段。
2. 建议在算法规格阶段确定时间步长、稳态容差、平均周期数和二维网格加密规则。依据是这些量决定非线性最优值的可复现性；风险是只写“搜索得到最优”无法支撑最大值结论；下一检查是补齐 `结果需求.md` 并生成 Q2/v2 同版本运行证据。

## Gate

- Can proceed to next stage: no
- Human decision needed: 审阅论文正文的保留/省略边界；补齐写作得分点预测后重新运行 `paper_model_distillation_check.py --for-paper`。
