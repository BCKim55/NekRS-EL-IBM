// NekRS-EL-IBM: one settings slide + one result slide per analysis
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");
const SKILL = process.env.PPTX_SKILL; // Claude pptx skill dir (apply_theme.js)
const { applyTheme } = require(path.join(SKILL, "scripts/apply_theme.js"));
const REPO = "/Users/byeong-cheonkim/Documents/NekRS-EL-IBM";
const IMG = process.env.DECK_IMG || path.join(__dirname, "deck_img"); // cropped panels (see docs/NekRS-EL-IBM_analyses.pptx)
const OUT = path.join(REPO, "docs", "NekRS-EL-IBM_analyses.pptx");

const THEME = {
  name: "NekRS-EL-IBM",
  headFontFace: "Calibri",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1C2B33", lt1: "FFFFFF", dk2: "0F4C5C", lt2: "EEF4F5",
    accent1: "0F4C5C", accent2: "E36414", accent3: "5B8E7D", accent4: "9A031E", accent5: "FB8B24", accent6: "4A6FA5",
    hlink: "0F4C5C", folHlink: "5B8E7D",
  },
};
const HEX = THEME.colors;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "NekRS-EL-IBM 해석 설정과 결과";
pres.author = "BC Kim";
const C = pres.SchemeColor;

// ---------- layouts ----------
pres.defineSlideMaster({
  title: "TITLE", background: { color: HEX.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.3, w: 11.7, h: 1.4, fontSize: 40, bold: true, color: C.background1, margin: 0 } } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.8, w: 11.7, h: 1.6, fontSize: 18, color: C.background2, margin: 0 } } },
  ],
});
pres.defineSlideMaster({
  title: "SECTION", background: { color: HEX.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.9, w: 11.7, h: 1.1, fontSize: 36, bold: true, color: C.background1, margin: 0 } } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.1, w: 11.7, h: 1.4, fontSize: 16, color: C.background2, margin: 0 } } },
  ],
});
pres.defineSlideMaster({
  title: "CONTENT", background: { color: HEX.lt1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 12.3, h: 0.8, fontSize: 26, bold: true, color: C.text2, margin: 0, valign: "middle" } } },
    { text: { text: "NekRS-EL-IBM  ·  해석 설정과 결과", options: { x: 0.5, y: 7.05, w: 6, h: 0.3, fontSize: 9, color: "8A9AA3", margin: 0 } } },
  ],
  slideNumber: { x: 12.3, y: 7.05, w: 0.6, h: 0.3, fontSize: 9, color: "8A9AA3" },
});

// ---------- helpers ----------
let section = "";
function content(title) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: section });
  s.addText(title, { placeholder: "title" });
  return s;
}
function kv(slide, rows, x, y, w, opts = {}) {
  const fs_ = opts.fontSize || 11.5;
  const c0 = opts.colW || [1.7, w - 1.7];
  slide.addTable(rows.map(([k, v]) => [
    { text: k, options: { bold: true, color: C.text2, fill: { color: C.background2 }, fontSize: fs_, valign: "middle" } },
    { text: v, options: { color: C.text1, fontSize: fs_, valign: "middle" } },
  ]), { x, y, w, colW: c0, border: { type: "solid", pt: 0.5, color: "C9D6DA" }, margin: 0.05, rowH: opts.rowH || 0.36, autoPage: false });
}
function grid(slide, header, rows, x, y, w, colW, opts = {}) {
  const fs_ = opts.fontSize || 11;
  const h = header.map((t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 }, fontSize: fs_, align: "center", valign: "middle" } }));
  const body = rows.map((r, i) => r.map((t, j) => ({ text: String(t), options: { fontSize: fs_, color: C.text1, align: j === 0 ? "left" : "center", valign: "middle", bold: j === 0, fill: { color: i % 2 ? C.background2 : C.background1 } } })));
  slide.addTable([h, ...body], { x, y, w, colW, border: { type: "solid", pt: 0.5, color: "C9D6DA" }, margin: 0.05, rowH: opts.rowH || 0.34, autoPage: false });
}
function stat(slide, x, y, w, value, label, color = C.accent2) {
  slide.addShape(pres.ShapeType.roundRect, { x, y, w, h: 1.25, fill: { color: C.background2 }, line: { color: C.background2 }, rectRadius: 0.08 });
  slide.addText(value, { x: x + 0.1, y: y + 0.08, w: w - 0.2, h: 0.65, fontSize: 24, bold: true, color, align: "center", margin: 0, isTextBox: true, fit: "shrink" });
  slide.addText(label, { x: x + 0.1, y: y + 0.72, w: w - 0.2, h: 0.48, fontSize: 10.5, color: C.text1, align: "center", valign: "top", margin: 0, isTextBox: true });
}
function bullets(slide, items, x, y, w, h, fontSize = 13) {
  slide.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 14 }, breakLine: i < items.length - 1, paraSpaceAfter: 6 } })),
    { x, y, w, h, fontSize, color: C.text1, valign: "top", margin: 0, isTextBox: true });
}
function caption(slide, text, x, y, w) {
  slide.addText(text, { x, y, w, h: 0.3, fontSize: 9.5, italic: true, color: "6B7C85", margin: 0, isTextBox: true });
}
function img(slide, file, x, y, w, h) {
  const p = fs.existsSync(file) ? file : path.join(IMG, file);
  slide.addImage({ path: p, x, y, w, h, objectName: path.basename(file) });
}
function pickColor(k) { return k; }

// Schematic: a 2D plan of the domain with labelled faces, bodies, particle clouds and arrows.
// spec: { aspect (w/h), faces: {top,bottom,left,right}, bodies:[{cx,cy,r,color,label}], clouds:[{cx,cy,w,h,n,color,r}],
//         arrows:[{x1,y1,x2,y2,label,color}], band:{x,y,w,h,label}, dashed:bool, caption, shape:'rect'|'circle' }
// cx, cy, w, h are fractions of the domain (origin bottom-left).
function schematic(slide, spec, area = { x: 0.5, y: 1.3, w: 5.9, h: 5.3 }) {
  const asp = spec.aspect || 1;
  const pad = 0.75;
  let bw = area.w - 2 * pad, bh = bw / asp;
  if (bh > area.h - 2 * pad - 0.4) { bh = area.h - 2 * pad - 0.4; bw = bh * asp; }
  const bx = area.x + (area.w - bw) / 2, by = area.y + pad + (area.h - 2 * pad - 0.4 - bh) / 2;
  const X = (f) => bx + f * bw, Y = (f) => by + (1 - f) * bh;
  const shape = spec.shape === "circle" ? pres.ShapeType.ellipse : pres.ShapeType.rect;
  slide.addShape(shape, { x: bx, y: by, w: bw, h: bh, fill: { color: "F6FAFB" }, line: { color: HEX.dk2, width: 1.5, dashType: spec.dashed ? "dash" : "solid" } });
  if (spec.band) {
    const b = spec.band;
    slide.addShape(pres.ShapeType.rect, { x: X(b.x), y: Y(b.y + b.h), w: b.w * bw, h: b.h * bh, fill: { color: "D9E8EB" }, line: { color: HEX.accent3, width: 0.75, dashType: "dash" } });
    if (b.label) slide.addText(b.label, { x: X(b.x), y: Y(b.y + b.h) - 0.28, w: b.w * bw, h: 0.26, fontSize: 9, color: HEX.accent3, align: "center", margin: 0, isTextBox: true });
  }
  const f = spec.faces || {};
  const fl = { fontSize: 10, color: HEX.dk2, margin: 0, isTextBox: true };
  if (f.top) slide.addText(f.top, { x: bx, y: by - 0.32, w: bw, h: 0.28, align: "center", ...fl });
  if (f.bottom) slide.addText(f.bottom, { x: bx, y: by + bh + 0.04, w: bw, h: 0.28, align: "center", ...fl });
  if (f.left) slide.addText(f.left, { x: bx - pad - 0.05, y: by, w: pad, h: bh, align: "center", valign: "middle", ...fl });
  if (f.right) slide.addText(f.right, { x: bx + bw + 0.05, y: by, w: pad, h: bh, align: "center", valign: "middle", ...fl });
  for (const c of spec.clouds || []) {
    const n = c.n || 40, r = c.r || 0.045;
    let seed = 7;
    const rnd = () => { seed = (seed * 9301 + 49297) % 233280; return seed / 233280; };
    for (let i = 0; i < n; i++) {
      const px = X(c.cx - c.w / 2 + rnd() * c.w), py = Y(c.cy - c.h / 2 + rnd() * c.h);
      slide.addShape(pres.ShapeType.ellipse, { x: px - r / 2, y: py - r / 2, w: r, h: r, fill: { color: c.color || HEX.accent2 }, line: { color: c.color || HEX.accent2, width: 0 } });
    }
    if (c.label) slide.addText(c.label, { x: X(c.cx - c.w / 2) - 0.3, y: Y(c.cy + c.h / 2) - 0.3, w: c.w * bw + 0.6, h: 0.26, fontSize: 9.5, color: HEX.accent2, align: "center", margin: 0, isTextBox: true });
  }
  for (const b of spec.bodies || []) {
    const d = 2 * b.r * bw;
    slide.addShape(pres.ShapeType.ellipse, { x: X(b.cx) - d / 2, y: Y(b.cy) - d / 2, w: d, h: d, fill: { color: b.color || "B5BEC4" }, line: { color: HEX.dk1, width: 1 } });
    if (b.label) slide.addText(b.label, { x: X(b.cx) - 1.0, y: Y(b.cy) + d / 2 + 0.03, w: 2.0, h: 0.26, fontSize: 9.5, color: HEX.dk1, align: "center", margin: 0, isTextBox: true });
  }
  for (const a of spec.arrows || []) {
    const x1 = X(a.x1), y1 = Y(a.y1), x2 = X(a.x2), y2 = Y(a.y2);
    slide.addShape(pres.ShapeType.line, { x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1) || 0.01, h: Math.abs(y2 - y1) || 0.01,
      flipH: x2 < x1, flipV: y2 < y1, line: { color: a.color || HEX.accent4, width: 1.75, endArrowType: "triangle" } });
    const lw = a.lw || 1.8;
    if (a.label) slide.addText(a.label, { x: (x1 + x2) / 2 - lw / 2, y: Math.min(y1, y2) - 0.3 + (a.dy || 0), w: lw, h: 0.26, fontSize: 9.5, color: a.color || HEX.accent4, align: "center", margin: 0, isTextBox: true });
  }
  if (spec.caption) caption(slide, spec.caption, area.x, area.y + area.h - 0.3, area.w);
}

function setupSlide(title, spec, rows, opts = {}) {
  const s = content(title);
  schematic(s, spec);
  kv(s, rows, 6.8, 1.3, 6.0, { rowH: opts.rowH, fontSize: opts.fontSize });
  if (opts.note) s.addText(opts.note, { x: 6.8, y: 6.25, w: 6.0, h: 0.6, fontSize: 10.5, italic: true, color: "6B7C85", margin: 0, isTextBox: true, valign: "top" });
  return s;
}

// ---------- title ----------
section = "개요";
pres.addSection({ title: section });
{
  const s = pres.addSlide({ masterName: "TITLE", sectionTitle: section });
  s.addText("NekRS-EL-IBM 해석 설정과 결과", { placeholder: "title" });
  s.addText("nekRS + ppiclF 오일러–라그랑주(EL) · 경계 매립법(IBM) 결합 코드\n검증 V1–V15와 Pilobolus 포자낭 발사 해석 — 해석별 설정 1장, 결과 1장\n2026-10-28", { placeholder: "body" });
}
{
  const s = content("진행한 해석 목록");
  const rows = [
    ["V1", "트레이서 강체회전 (일방향 결합)", "검증", "반올림 수준 (2e-14)"],
    ["V2", "두 방향 결합 보존 수지", "검증", "중력 충격량의 1e-3"],
    ["V3", "디바이스 입자 경로 vs ppiclF 경로", "검증", "0.2–0.6 %"],
    ["V4", "분무 파셀 2D 슬롯 제트, 두 경로", "검증", "7e-8; 0.8 %"],
    ["V5", "EL–IBM 접촉과 운동량 교환", "검증", "운동량 수지 0.3 %"],
    ["V6", "ppiclF 패치 전·후 회귀", "검증", "패치판 1–3 랭크 정확"],
    ["V7", "nekRS+ppiclF vs Nek5000+ppiclF", "교차 검증", "0.2–1.2 %"],
    ["V8", "TAB 변형·분열", "검증+실험 비교", "주기 0.5 %; Pilch & Erdman"],
    ["V9", "파셀 평면 제트 vs Bickley", "실험·이론 비교", "u_c 3.5 % 이내"],
    ["V10", "고정 구 Re 10 (NekIBM simple)", "실험·이론 비교", "4.09 / 4.69 vs 4.15–4.30"],
    ["V11", "고정 구 Re 100 vs Johnson & Patel", "실험·이론 비교", "C_d +18 % → +10 %"],
    ["V12", "침강 유리구 vs Mordant & Pinton", "실험 비교", "u_t −6 % → −2.5 %"],
    ["V13", "EL·IBM 필터 폭 분리", "검증", "동일 폭: 결과 재현"],
    ["V14", "처방 운동 구 vs 고정 구", "검증", "C_d 1.257 vs 1.289"],
    ["V15", "부착 접촉 (sticking)", "검증", "107/200 포획, 수지 동일"],
    ["P1", "Pilobolus 단상 발사 (처방·자유비행)", "응용", "공기 감속 111–125 m/s²"],
    ["P2", "Pilobolus 수액 액적 부착", "응용", "부착 질량 90 %, 항력 −36 %"],
  ];
  grid(s, ["#", "해석", "종류", "핵심 결과"], rows, 0.5, 1.25, 12.3, [0.7, 5.3, 2.0, 4.3], { fontSize: 10.5, rowH: 0.31 });
}

// ---------- V1–V9: EL coupling ----------
section = "1. EL 결합 검증 (V1–V9)";
pres.addSection({ title: section });
{
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: section });
  s.addText("1. EL 결합 검증 (V1–V9)", { placeholder: "title" });
  s.addText("nekRS ↔ ppiclF 배관, 두 방향 결합 보존, 디바이스 경로, 분무 파셀, 충돌, 교차 코드, TAB, Bickley 제트", { placeholder: "body" });
}

// V1
setupSlide("V1 설정 — 트레이서 강체회전 (일방향 결합)", {
  shape: "circle", aspect: 1, faces: { top: "ethier.re2: 단위 구, 32 요소, N = 7" },
  clouds: [{ cx: 0.5, cy: 0.5, w: 0.7, h: 0.7, n: 60, r: 0.05, color: HEX.accent2, label: "" }],
  arrows: [{ x1: 0.82, y1: 0.35, x2: 0.82, y2: 0.65, label: "u = (−ωy, ωx), ω = 1", dy: -0.1 }],
  caption: "유동을 강체회전으로 고정하고 트레이서 1000개를 한 바퀴 돌림 (1·2·4·8 랭크)",
}, [
  ["격자", "nekRS ethier.re2 (단위 구, 32 요소), N = 7"],
  ["유동", "강체회전 u = −ωy, v = ωx (ω = 1)로 고정"],
  ["입자", "트레이서 1000개, 무작위 위치, 일방향"],
  ["시간", "dt = 2π/600, 600 스텝 (1회전)"],
  ["랭크", "1, 2, 4, 8"],
  ["기준해", "속도가 선형이라 GLL 보간은 정확; ppiclF RK3는 이 RHS에서 전진 오일러와 같아 이산 해석해 존재"],
  ["검사 항목", "보간, 요소 소유, 랭크 간 입자 이동"],
], { rowH: 0.5 });
{
  const s = content("V1 결과 — 반올림 수준에서 정확, 랭크 수 무관");
  stat(s, 0.5, 1.4, 3.9, "1.9e-14", "이산 해석해와의 최대 거리 (step 550, 8 랭크)");
  stat(s, 4.7, 1.4, 3.9, "2.1e-14", "같은 값 (step 600); 1·2·4 랭크 동일");
  stat(s, 8.9, 1.4, 3.9, "0", "손실된 트레이서");
  bullets(s, [
    "errExact (정확한 회전과의 차) 1.7e-2는 ppiclF 자체의 O(dt) 오일러 오차로 1회전 예측치와 일치 — 결합의 오차가 아님",
    "보간·요소 소유·이동이 랭크 수에 관계없이 반올림 수준에서 정확",
    "힘 투영이나 두 방향 결합 항은 시험하지 않음 (→ V2)",
  ], 0.5, 3.1, 12.3, 2.2);
  s.addShape(pres.ShapeType.roundRect, { x: 0.5, y: 5.4, w: 12.3, h: 1.3, fill: { color: C.background2 }, line: { color: C.background2 }, rectRadius: 0.08 });
  s.addText([{ text: "의미  ", options: { bold: true, color: C.text2 } }, { text: "nekRS 요소 배치를 ppiclF 겹침 격자로 넘기는 배관과 입자 재분배가 올바르게 동작한다. 이후 모든 두 방향 검증의 바탕." }],
    { x: 0.7, y: 5.5, w: 11.9, h: 1.1, fontSize: 13, color: C.text1, margin: 0, isTextBox: true, valign: "middle" });
}

// V2
setupSlide("V2 설정 — 두 방향 결합 보존 수지 (settling_twoway)", {
  aspect: 1, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  clouds: [{ cx: 0.5, cy: 0.5, w: 0.9, h: 0.9, n: 120, r: 0.035, color: HEX.accent2 }],
  arrows: [{ x1: 0.08, y1: 0.9, x2: 0.08, y2: 0.65, label: "g = −1 (z)", dy: 0.02 }],
  caption: "완전 주기 단위 상자, 정지 입자 10 000개가 중력으로 침강 (6³ 요소, N = 5)",
}, [
  ["영역·격자", "완전 주기 단위 상자, 6 × 6 × 6 요소, N = 5"],
  ["유체", "ρ = 1, ν = 0.01"],
  ["입자", "점입자 10 000개, d_p 0.01, ρ_p 100 (φ_p 0.52 %, 질량부하 0.52), 정지 상태에서 출발"],
  ["결합", "가우시안 투영 FWHM 0.15, 운동량 소스, QTL 발산 (qtl = 1)"],
  ["시간", "dt = 1e-3, 200 스텝 (t = 0.2)"],
  ["랭크", "1, 4"],
  ["검사 수지", "proj (투영 보존), fluid (소스 수령), total (전체 운동량), div / divErr (부과·해석 발산)"],
], { rowH: 0.5 });
{
  const s = content("V2 결과 — 세 가지 전달이 각각 보존 수지를 닫음");
  grid(s, ["수지 (t = 0.2)", "값", "의미"], [
    ["proj: ∫F dV − Σf_p", "3.6e-4", "가우시안 절단(1e-4)과 GLL 구적의 보존 오차"],
    ["fluid: P_f − P_f0 − ∫∫FU", "1.6e-3 (중력 충격량 대비)", "유체 솔버가 소스를 그대로 받음"],
    ["total: P_f + P_p − P_0 − Σm g t", "−1.0e-3", "외력은 중력뿐; 명시적 결합 지연만큼 손실"],
    ["div, divErr", "2.0e-4, 6.0e-5", "부과 발산과 해석 발산의 차 = nekRS 자체 오차"],
    ["1 랭크 vs 4 랭크", "출력 자릿수까지 동일", "병렬 독립"],
  ], 0.5, 1.3, 12.3, [3.6, 3.0, 5.7], { fontSize: 11.5, rowH: 0.45 });
  stat(s, 0.5, 4.3, 3.9, "1e-3", "전체 운동량 손실 / 중력 충격량");
  stat(s, 4.7, 4.3, 3.9, "3.6e-4", "투영 보존 오차 (절단 한계)");
  stat(s, 8.9, 4.3, 3.9, "동일", "1 vs 4 랭크");
  bullets(s, ["항력 법칙이나 필터 폭이 물리적으로 맞는지는 보여주지 않음 (→ V7, V9)"], 0.5, 5.9, 12.3, 0.8, 12);
}

// V3
setupSlide("V3 설정 — 디바이스 입자 경로 vs ppiclF 경로 (settling_gpu)", {
  aspect: 1, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  clouds: [{ cx: 0.5, cy: 0.5, w: 0.9, h: 0.9, n: 120, r: 0.035, color: HEX.accent6 }],
  arrows: [{ x1: 0.08, y1: 0.9, x2: 0.08, y2: 0.65, label: "g = −1 (z)", dy: 0.02 }],
  caption: "V2와 동일한 문제를 src/elgpu.hpp (nekRS lpm_t, 자체 커널)로 해석",
}, [
  ["문제", "V2와 동일 (상자, 입자, 필터 폭, qtl = 1, 200 스텝, 1·4 랭크)"],
  ["ppiclF 경로", "ppiclF 비닝·고스트·가우시안 합 투영, RK3"],
  ["디바이스 경로", "nekRS lpm_t + 자체 커널: L2 투영 + 확산 필터, QTL, SSP-RK3"],
  ["공유 코드", "입자 쪽은 없음 (독립 구현 두 개)"],
  ["백엔드", "OCCA serial (CPU); 커널은 @atomic만 사용"],
  ["검사 항목", "P_f, P_p, Σf, 발산, 투영 보존, 스텝 시간"],
], { rowH: 0.52 });
{
  const s = content("V3 결과 — 두 독립 구현이 1 % 이내에서 일치");
  grid(s, ["t = 0.2", "ppiclF 경로 (V2)", "디바이스 경로", "차이"], [
    ["유체 운동량 P_f", "−5.631e-2", "−5.666e-2", "0.6 %"],
    ["입자 운동량 P_p", "−4.851e-2", "−4.865e-2", "0.3 %"],
    ["항력 합 Σf", "−0.3410", "−0.3417", "0.2 %"],
    ["부과 발산 RMS", "2.01e-4", "1.98e-4", ""],
    ["투영 보존", "4.6e-4 (절단)", "3e-16 (반올림)", ""],
    ["전체 운동량 드리프트", "−1.0e-3", "−1.2e-3", ""],
    ["스텝 시간 1 / 4 랭크", "2.0 s / 0.9 s", "0.24 s / 0.09 s", "8–10x"],
  ], 0.5, 1.3, 7.6, [2.5, 1.9, 1.9, 1.3], { fontSize: 11, rowH: 0.42 });
  stat(s, 8.5, 1.3, 4.3, "0.2–0.6 %", "운동량·항력 차이");
  stat(s, 8.5, 2.75, 4.3, "8–10x", "디바이스 경로 속도 (CPU 직렬 기준)");
  bullets(s, [
    "남은 차이는 필터(가우시안 합 vs 확산)와 힘을 취하는 RK 단계",
    "GPU 백엔드에서는 아직 실행하지 않음 (OCCA serial만)",
  ], 8.5, 4.3, 4.3, 2.2, 12);
}

// V4
setupSlide("V4 설정 — 분무 파셀 2D 슬롯 제트 (spray_jet2d, _gpu)", {
  aspect: 2, faces: { top: "free-slip", bottom: "free-slip", left: "주기 (x)", right: "주기 (x)" },
  clouds: [{ cx: 0.3, cy: 0.5, w: 0.35, h: 0.2, n: 70, r: 0.03, color: HEX.accent2, label: "파셀 (2.4 액적/파셀)" }],
  arrows: [{ x1: 0.05, y1: 0.5, x2: 0.14, y2: 0.5, label: "슬롯 0.02, u_jet = 1", dy: 0.35 }],
  caption: "2 × 1 × 0.1 상자 (40 × 20 × 1 요소, N = 5), z 주기 — 준 2D",
}, [
  ["영역·격자", "2 × 1 × 0.1, 40 × 20 × 1 요소, N = 5; x·z 주기, y free-slip"],
  ["기체", "ρ = 1, ν = 5e-4"],
  ["액체", "ρ 1000, d_p 0.002 (τ_p 0.44, Re_p0 4), 슬롯 폭 0.02 (x = 0.1), u_jet 1, 슬롯 내 φ 0.01 (질량부하 10)"],
  ["주입", "스텝당 파셀 10개 × 2.4 액적, 두 경로가 같은 결정론적 생성기 사용"],
  ["시간", "dt = 2e-3, 1000 스텝 (t = 2), 중력 없음, 4 랭크"],
  ["기준", "(i) 단일 액적 vs 항력 ODE (Python), (ii) 주입 운동량 (정확히 알려짐), (iii) 상대 경로"],
], { rowH: 0.55 });
{
  const s = content("V4 결과 — 두 경로 모두 액적 ODE·운동량 수지를 만족");
  img(s, path.join(REPO, "docs/img/spray_jet2d.png"), 0.5, 1.25, 3.9, 5.37);
  grid(s, ["", "ppiclF 경로", "디바이스 경로"], [
    ["단일 액적 vs 항력 ODE, 최대 오차/U0", "7e-8", "7e-8"],
    ["운동량 수지 t = 2 (최대)", "+8.6e-4 (2.8e-3)", "+9.4e-4 (2.8e-3)"],
    ["유체 수지 P_gas − ∫FU", "1.6e-4", "1.6e-4"],
    ["P_gas, P_parcels (t = 2)", "1.072e-2, 2.931e-2", "1.081e-2, 2.923e-2"],
    ["중심선 u_x 최대", "0.8333", "0.8336"],
    ["중심선 차 (max / rms)", "—", "0.0055 / 0.0013"],
    ["스텝 시간 (1만 파셀)", "1.25 s", "0.19 s"],
  ], 4.7, 1.3, 8.1, [3.7, 2.2, 2.2], { fontSize: 11, rowH: 0.4 });
  bullets(s, [
    "액적 방정식은 RK3 정확도로 적분, 운동량은 주입량의 1e-3 이내 보존, 두 경로 차 0.8 % (운동량) / 0.7 % (중심선)",
    "dt = 5e-3 실행은 유도 기체 속도(0.83 u_jet)가 nekRS 이류 CFL을 넘어 t = 0.5에서 발산 — 유체 솔버 한계, 결합 문제 아님",
  ], 4.7, 5.3, 8.1, 1.5, 11);
  caption(s, "위: 두 경로의 기체 중심선 속도와 파셀 분포 (docs/img/spray_jet2d.png)", 0.5, 6.65, 4.0);
}

// V5
setupSlide("V5 설정 — EL–IBM 접촉과 운동량 교환 (el_ibm_collision)", {
  aspect: 1, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  bodies: [{ cx: 0.5, cy: 0.5, r: 0.15, label: "IBM 구 R 0.15, ρ_b/ρ_f 4 (자유)" }],
  clouds: [{ cx: 0.18, cy: 0.5, w: 0.2, h: 0.2, n: 60, r: 0.03, color: HEX.accent2, label: "EL 200개, d_p 0.01" }],
  arrows: [{ x1: 0.3, y1: 0.72, x2: 0.42, y2: 0.72, label: "u = 0.5", dy: 0.02 }],
  caption: "완전 주기 단위 상자, 8³ 요소, N = 5 — 정지한 자유 IBM 구에 입자 구름을 발사",
}, [
  ["영역·격자", "완전 주기 단위 상자, 8³ 요소, N = 5 (GLL ≈ 0.025)"],
  ["유체", "ρ = 1, μ = 0.01, dt 5e-3, 240 스텝 (t = 1.2), 4 랭크"],
  ["IBM 구", "R 0.15 (d = 0.3, 격자 12칸), ρ_b/ρ_f 4, 마커 간격 0.025, 가우시안 FWHM 0.06"],
  ["EL 입자", "200개, d_p 0.01 (격자 0.4칸), ρ_p 2000, +x로 0.5"],
  ["접촉", "연성 구 k_n 500, 반발계수 0.9, 입자 서브스텝 10 (접촉시간 4.5e-3 vs dt 5e-3)"],
  ["기준", "외력 없음 → P_fluid + P_EL + P_body = P_0; gap/d_p, 접촉 수, 반발 수"],
], { rowH: 0.55 });
{
  const s = content("V5 결과 — 입자가 구 안으로 들어가지 않고 운동량이 0.3 % 이내로 닫힘");
  img(s, "collision_snaps.png", 0.5, 1.25, 7.4, 0.87);
  img(s, "collision_vel.png", 0.5, 2.3, 3.9, 2.53);
  img(s, path.join(REPO, "cases/el_ibm_collision/snapshots/frames/frame_0030.png"), 4.6, 2.3, 3.4, 2.55);
  grid(s, ["t = 1.2", "값"], [
    ["구 속도", "0 → 0.271"],
    ["P_f / P_EL / P_b", "3.74e-2 / 5.55e-2 / 1.15e-2"],
    ["(ΣP − P_0)/P_0", "−3.2e-3"],
    ["최소 gap/d_p", "0.109 (내부 입자 없음)"],
    ["반발 입자", "4 / 200"],
    ["마커 슬립 RMS", "4.8e-3 (구 속도의 2 %)"],
  ], 8.3, 1.3, 4.5, [2.0, 2.5], { fontSize: 11, rowH: 0.4 });
  bullets(s, [
    "손실 0.3 %는 강성 접촉에서의 명시적 결합 지연; 마커 간격과 무관",
    "구가 원방장을 크게 바꾸지 않아 교환의 검증이지 항력 검증은 아님",
  ], 8.3, 4.2, 4.5, 2.2, 11.5);
  caption(s, "위: 자오면 스냅샷(t = 0, 0.4, 0.8, 1.2). 아래 왼쪽: 구·EL 속도. 오른쪽: ParaView (φ_p 단면, t = 0.6)", 0.5, 5.1, 7.5);
}

// V6
setupSlide("V6 설정 — ppiclF 라이브러리 패치 전·후 회귀 (ppiclf_regression.py)", {
  aspect: 1.3, faces: { top: "ppiclF examples/stokes_2d (랭크 1, 2, 3)" },
  clouds: [{ cx: 0.5, cy: 0.6, w: 0.9, h: 0.5, n: 90, r: 0.035, color: HEX.accent2, label: "랭크당 250개" }],
  arrows: [{ x1: 0.08, y1: 0.45, x2: 0.08, y2: 0.15, label: "g = 9.8, Stokes τ = 1/9.8", dy: 0.02 }],
  caption: "정확해 v(t) = −gτ(1 − e^(−t/τ))를 입자마다 VTU 태그로 대조",
}, [
  ["대상", "patches/ppiclF/0001–0005 (이동 후 RK3 경계, 겹침 격자 노드 검사, 투영 종료, 주기 고스트 중복 제거, 필터-주기 abort)"],
  ["빌드", "원본·패치판을 같은 방식으로 (gfortran, -DPPICLC 없음)"],
  ["문제", "ppiclF 자체 예제 stokes_2d: 랭크당 250개, Stokes 항력 τ_p 1/9.8, g 9.8, 정지 출발, dt 1e-4, 10 000 스텝"],
  ["랭크", "1, 2, 3; ppiclF tests/test1 (보간·투영 진단)도 비교"],
  ["기준", "정확해와 그 적분, 입자별 (VTU TIME은 스텝 시작 시각, 상태는 dt 뒤)"],
], { rowH: 0.62 });
{
  const s = content("V6 결과 — 원본은 다중 랭크에서 RK 단계를 잃고, 패치판은 정확");
  grid(s, ["라이브러리", "랭크", "v 최대 오차", "변위 최대 오차", "1e-5 이상 벗어난 입자"], [
    ["원본", "1", "3.7e-8", "3.4e-8", "0"],
    ["원본", "2, 3", "1.8e-3", "1.7e-3", "1 / 500, 1 / 750"],
    ["패치판", "1, 2, 3", "3.7e-8", "3.4e-8", "0"],
  ], 0.5, 1.3, 12.3, [2.4, 1.6, 2.4, 2.6, 3.3], { fontSize: 12, rowH: 0.45 });
  stat(s, 0.5, 3.4, 3.9, "3.7e-8", "패치판 속도 오차, 모든 랭크 수");
  stat(s, 4.7, 3.4, 3.9, "1 / 500", "원본에서 이동 후 RK 단계를 잃는 입자 (패치 0001)");
  stat(s, 8.9, 3.4, 3.9, "동일", "test1 진단 (TimeI 9.3e-9, Prjct 1.3e-8)");
  bullets(s, [
    "패치 0001이 고치는 결함은 2 랭크 이상에서만 나타남; 나머지 패치는 조건이 아니면 결과를 바꾸지 않음",
    "이 회귀로 V7의 Nek5000 비교에서는 원본 라이브러리를 4 랭크로 써도 되는지(정확한 랭크 수) 판단",
  ], 0.5, 5.0, 12.3, 1.6, 12);
}

// V7
setupSlide("V7 설정 — nekRS + ppiclF vs Nek5000 + 원본 ppiclF (교차 코드)", {
  aspect: 1, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  clouds: [{ cx: 0.5, cy: 0.5, w: 0.9, h: 0.9, n: 120, r: 0.035, color: HEX.accent3 }],
  arrows: [{ x1: 0.08, y1: 0.9, x2: 0.08, y2: 0.65, label: "g = −1", dy: 0.02 }],
  caption: "V2의 상자를 ppiclF가 원래 쓰이던 Nek5000 v19 예제 배선으로 그대로 해석",
}, [
  ["기준 코드", "Nek5000 v19 + 원본(비패치) ppiclF, Nek5000 예제 배선 (GetProFldIJKEF 소스 F/(1−φ_p)/ρ, userchk에서 보간·IntegrateParticle), BDF2"],
  ["이 코드", "nekRS + ppiclF, BDF3/EXT, design.md §1의 소스 처리"],
  ["공통", "V2의 격자·dt·필터·입자 수 (무작위 위치만 다름), 4 랭크씩"],
  ["비교량", "P_f (유체 운동량), P_p (입자 운동량), Σf (항력 합) at t = 0.04, 0.12, 0.20"],
  ["자료", "cases/settling_twoway/nek5000_reference (입력, README, 로그)"],
], { rowH: 0.62 });
{
  const s = content("V7 결과 — t = 0.2에서 두 코드가 0.2–1.2 % 이내로 일치");
  img(s, path.join(REPO, "docs/img/validation.png"), 0.5, 1.25, 12.3, 3.53);
  grid(s, ["t", "코드", "P_f", "P_p", "Σf"], [
    ["0.04", "Nek5000 / nekRS", "−5.682e-3 / −5.333e-3", "−1.583e-2 / −1.552e-2", "−0.2306 / −0.2283"],
    ["0.12", "Nek5000 / nekRS", "−2.974e-2 / −2.919e-2", "−3.389e-2 / −3.362e-2", "−0.3297 / −0.3300"],
    ["0.20", "Nek5000 / nekRS", "−5.696e-2 / −5.631e-2", "−4.883e-2 / −4.851e-2", "−0.3404 / −0.3410"],
  ], 0.5, 4.95, 8.3, [0.7, 1.8, 2.0, 2.0, 1.8], { fontSize: 10.5, rowH: 0.36 });
  bullets(s, [
    "t = 0.2: P_f 1.2 %, P_p 0.7 %, Σf 0.2 %",
    "초기(t = 0.04) 6 % 차이는 시간적분(BDF2 vs BDF3/EXT)과 다른 무작위 위치",
    "운동량 수지: Nek5000 1 %, nekRS 0.1 %; 디바이스 경로(V3)도 같은 띠 안",
  ], 9.1, 4.95, 3.7, 2.0, 11);
}

// V8
setupSlide("V8 설정 — TAB 변형·분열 (tab_breakup)", {
  aspect: 1.6, faces: { top: "정지 기체 (ρ 1, ν 5e-5), 일방향" },
  bodies: [{ cx: 0.25, cy: 0.5, r: 0.07, color: HEX.accent2, label: "We 10 … 300" }],
  arrows: [{ x1: 0.4, y1: 0.5, x2: 0.62, y2: 0.5, label: "상대속도 U → 변형 y(t)", dy: 0.02 }, { x1: 0.68, y1: 0.5, x2: 0.9, y2: 0.5, label: "y = 1 → 분열", dy: 0.02, color: HEX.accent4 }],
  caption: "단일 액적 11개를 서로 다른 We로 정지 기체에 놓고 TAB 진동자 적분 (디바이스 경로)",
}, [
  ["액적", "d 0.002, σ 1e-4, ρ_l 1000, μ_l 1e-3 (Oh 0.07)"],
  ["We", "10, 12, 15, 18, 25, 35, 50, 80, 120, 200, 300"],
  ["시간", "dt 2e-3, 600 스텝, 일방향 결합, 1 랭크"],
  ["기준 (i)", "처음 3 스텝의 y, y', v vs Python SSP-RK3"],
  ["기준 (ii)", "감쇠 진동 주기 2π/√(ω² − ζ²), ω² = 8σ/(ρ_l r³)"],
  ["기준 (iii, iv)", "TAB 임계 We = 12; Pilch & Erdman (1987) T_ini = 1.9 (We − 12)^(−1/4)"],
], { rowH: 0.5 });
{
  const s = content("V8 결과 — 진동자는 정확, 분열 개시는 저 We에서 실험 곡선 위");
  img(s, path.join(REPO, "docs/img/tab_breakup.png"), 0.5, 1.25, 8.0, 2.91);
  grid(s, ["검사", "결과"], [
    ["처음 3 스텝 vs Python SSP-RK3", "6자리 일치"],
    ["We 10 진동 주기", "0.224 / 0.222 / 0.224 vs 해석 0.2230 (0.5 %)"],
    ["We 10, 12, 15 분열", "없음 (감쇠가 임계를 12 위로)"],
    ["T (We 18 → 120)", "1.17, 1.06, 1.00, 1.00, 0.95, 0.93"],
    ["Pilch & Erdman T_ini", "1.21, 1.00, 0.87, 0.77, 0.66, 0.59"],
    ["자식 크기 d_32/d", "0.40 (We 18) → 0.10 (We 120)"],
  ], 0.5, 4.35, 8.0, [3.0, 5.0], { fontSize: 10, rowH: 0.3 });
  bullets(s, [
    "저 We에서는 측정 개시 곡선 위, We 120에서는 1.6x 위 (TAB 포화)",
    "미해결: We ≥ 200에서 자식이 몇 스텝 만에 재분열 → τ_p < dt → 명시적 항력 불안정; dt를 줄이거나 음해법 항력 필요",
  ], 8.8, 1.3, 4.0, 5.4, 11);
}

// V9
setupSlide("V9 설정 — 파셀 평면 제트 vs Bickley (jet_bickley)", {
  aspect: 1.25, faces: { top: "free-slip", bottom: "free-slip", left: "유입 (5 % 동류)", right: "유출" },
  clouds: [{ cx: 0.35, cy: 0.5, w: 0.4, h: 0.1, n: 80, r: 0.025, color: HEX.accent2, label: "파셀 d_p 0.0005" }],
  arrows: [{ x1: 0.05, y1: 0.5, x2: 0.14, y2: 0.5, label: "슬롯 0.04, u_jet 10", dy: 0.35 }, { x1: 0.6, y1: 0.5, x2: 0.9, y2: 0.5, label: "Bickley 제트 u_c ∝ x^(−1/3)", dy: 0.02, color: HEX.accent3 }],
  caption: "열린 슬래브, 40 × 32 × 1 요소, N = 5; 희박(질량부하 0.01)하게 주입해 기체만 남김",
}, [
  ["영역·격자", "유입(동류 5 %)·유출, 옆면 free-slip, 40 × 32 × 1 요소, N = 5, ν 5e-4"],
  ["주입", "슬롯 0.04 (x = 0.5), u_jet 10, d_p 0.0005 (τ_p 0.016, 정지 거리 0.16), 슬롯 내 φ 5e-6 (질량부하 0.01)"],
  ["시간", "스텝당 파셀 10개, dt 2.5e-3, 8000 스텝 (t = 20, 정상), 마지막 11 700 파셀, 4 랭크"],
  ["대조군", "같은 J를 고정 가우시안 체적력으로 준 nekRS 단독 실행"],
  ["기준", "Bickley: u_c = 0.4543 (J²/νx)^(1/3), u = u_c sech²(…); J는 실행에서 측정한 값 (0.0173)"],
  ["민감도", "필터 폭 0.04 → 0.06, 스텝당 파셀 10 → 40"],
], { rowH: 0.55 });
{
  const s = content("V9 결과 — 파셀 제트가 Bickley의 감쇠와 sech² 분포를 재현");
  img(s, path.join(REPO, "docs/img/bickley.png"), 0.5, 1.25, 6.2, 4.51);
  grid(s, ["0.9 < x < 1.8", "파셀", "체적력"], [
    ["∫u(u − u_co)dy (창 내 일정)", "0.0169 (0.3 %)", "0.0173 (0.3 %)"],
    ["u_c 감쇠에서 구한 J", "0.0168", "0.0165"],
    ["가상 원점 x0", "0.572 (노즐 + 정지거리)", "0.325"],
    ["u_c vs Bickley 최대 오차", "3.5 %", "3.7 %"],
    ["반폭 vs Bickley 최대 오차", "13 %", "12 %"],
    ["x = 1.2 분포", "sech²", "sech²"],
  ], 7.0, 1.3, 5.8, [2.6, 1.7, 1.5], { fontSize: 10.5, rowH: 0.38 });
  bullets(s, [
    "파셀 제트는 체적력 제트보다 2.3 % 적게 운반 (0.5 %는 파셀에 남음)",
    "반폭 오차는 옆면 구속(두 실행 공유, 서로 6 % 차) — 결합 문제 아님",
    "민감도: J는 세 자리까지 동일, 감쇠 상수 1–2 %, 필터는 가상 원점만 이동 (0.57 → 0.51)",
  ], 7.0, 4.7, 5.8, 2.1, 10.5);
  caption(s, "scripts/jet_bickley.py: 중심선 속도, 반폭, x = 1.2 분포 (docs/img/bickley.png)", 0.5, 5.85, 6.2);
}

// ---------- V10–V12: IBM ----------
section = "2. IBM 검증 (V10–V12)";
pres.addSection({ title: section });
{
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: section });
  s.addText("2. IBM 검증 (V10–V12)", { placeholder: "title" });
  s.addText("NekIBM의 입력으로 고정 구 Re 10·Re 100, 침강 유리구 — 노트북 격자와 Tuolumne(MI300A) 고해상도", { placeholder: "body" });
}

// V10
setupSlide("V10 설정 — 고정 구 Re 10, NekIBM simple 상자 (ibm_fixed_sphere, _fine)", {
  aspect: 1, faces: { top: "유출", bottom: "유입 u = 1 (y)", left: "주기 (x, z)", right: "주기 (x, z)" },
  bodies: [{ cx: 0.5, cy: 0.5, r: 0.125, label: "구 d = 0.01, 마커 1963개" }],
  arrows: [{ x1: 0.5, y1: 0.08, x2: 0.5, y2: 0.3, label: "u_in", dy: 0.05 }],
  caption: "(4d)³ 상자, 20³ 요소, N = 5 (GLL ≈ d/25); 가우시안 폭 두 가지 (0.08 d, 0.04 d)",
}, [
  ["영역·격자", "NekIBM IBM_test/simple: (4d)³, 20³ 요소, N = 5 (GLL ≈ d/25)"],
  ["유체", "유입 u = 1, 유출 (y), x·z 주기, ν 1e-3 (Re = 10)"],
  ["IBM", "마커 1963개 (d/25); FWHM 0.08 d (GLL 2칸) vs _fine 0.04 d (1칸, 마커 간격 절반)"],
  ["시간", "dt 2.5e-4, 200 스텝 (U t/d = 5), 4 랭크; Tuolumne 4 MI300A 재실행"],
  ["측정", "C_d = F_y/(ρ u² π d²/8) (마커 힘 합), 횡력, 마커 슬립"],
  ["기준", "비구속 구 Re 10: Schiller–Naumann 4.15, 표준 항력 4.30; 상자는 비구속 아님 (측면 주기 4d, φ 0.8 %)"],
], { rowH: 0.55 });
{
  const s = content("V10 결과 — 두 가우시안 폭이 비구속 기준값을 양쪽에서 감쌈");
  img(s, "re10_curve.png", 0.5, 1.25, 4.8, 4.04);
  img(s, "re10_slices.png", 5.5, 1.25, 7.3, 4.02);
  grid(s, ["FWHM", "C_d (U t/d = 5)", "마커 슬립 / u_in", "횡력"], [
    ["0.08 d (GLL 2칸)", "4.694", "18 %", "1e-10"],
    ["0.04 d (GLL 1칸)", "4.090", "31 %", "1e-10"],
    ["Tuolumne, 4 MI300A, 0.08 d", "4.6937 (비트 동일)", "", ""],
  ], 0.5, 5.45, 7.4, [2.6, 1.9, 1.7, 1.2], { fontSize: 10.5, rowH: 0.33 });
  bullets(s, [
    "넓은 가우시안은 유효 몸체를 키우고(정체 영역이 원 밖까지), 한 칸 폭은 no-slip을 덜 강제(슬립 31 %)",
    "항력 수렴은 필터를 좁히는 게 아니라 격자를 세밀하게 해야 함 (→ V11); 대칭대로 횡력 소멸",
  ], 8.1, 5.45, 4.7, 1.5, 10.5);
}

// V11
setupSlide("V11 설정 — 고정 구 Re 100 vs Johnson & Patel (ibm_sphere_re100, hires)", {
  aspect: 0.6, faces: { top: "유출", bottom: "유입 u = 1 (y)", left: "주기 (x, z)", right: "주기 (x, z)" },
  bodies: [{ cx: 0.5, cy: 0.4, r: 0.0625, label: "구 d = 0.01" }],
  arrows: [{ x1: 0.5, y1: 0.05, x2: 0.5, y2: 0.22, label: "u_in", dy: 0.5 }],
  caption: "NekIBM fixsphere/Re100 상자 8d × 15d × 8d; 노트북 격자와 NekIBM 격자(hires)",
}, [
  ["영역", "8d × 15d × 8d, d 0.01, 유입·유출 (y), x·z 주기, ν 1e-4 (Re 100)"],
  ["노트북", "16 × 30 × 16 요소, N 5 (d/10), 마커 d/17, FWHM 0.12 d, dt 2e-4, 500 스텝, 4 랭크 (CPU)"],
  ["hires", "32 × 60 × 32 (NekIBM 격자), N 7 (d/28, 모서리 d/60), 마커 d/50 (7854개), FWHM 0.04 d, dt 5e-5, 2000 스텝, 8 MI300A (1.1 s/스텝, 41 분)"],
  ["종료", "둘 다 U t/d = 10; 마지막 스텝에 z = 0 평면 저장 (fieldDump.hpp)"],
  ["기준", "Johnson & Patel (1999): C_d 1.09, 정상 와 길이 0.87 d (정상 상태)"],
], { rowH: 0.65 });
{
  const s = content("V11 결과 — 해상도를 올리면 항력 오차가 반으로 (+18 % → +10 %)");
  img(s, "re100_curve.png", 0.5, 1.25, 4.9, 4.03);
  img(s, "re100_slices.png", 5.6, 1.25, 5.26, 4.03);
  grid(s, ["U t/d = 10", "노트북", "hires", "J & P"], [
    ["C_d", "1.286 (+18 %)", "1.201 (+10 %)", "1.09"],
    ["정상 와 길이", "1.00 d", "0.95 d", "0.87 d"],
    ["마커 슬립 / u_in", "3.9 %", "2.8 %", ""],
    ["dC_d/d(Ut/d) 끝", "", "−0.002", "0"],
  ], 0.5, 5.45, 7.0, [2.2, 1.7, 1.7, 1.4], { fontSize: 10.5, rowH: 0.33 });
  bullets(s, [
    "가우시안 경계 폭이 주 오차; 세밀화하면 Johnson & Patel 쪽으로 수렴",
    "남은 10 %의 일부는 U t/d = 10에서 아직 발달 중인 후류, 일부는 경계 폭; 8d 주기(φ 0.5 %)의 작은 폐색",
  ], 7.7, 5.45, 5.1, 1.5, 10.5);
}

// V12
setupSlide("V12 설정 — 침강 유리구 vs Mordant & Pinton (ibm_settling_sphere)", {
  aspect: 0.4, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  bodies: [{ cx: 0.5, cy: 0.8, r: 0.07, label: "d = 1/6, ρ_p/ρ_f 2.56" }],
  arrows: [{ x1: 0.8, y1: 0.7, x2: 0.8, y2: 0.45, label: "g = 9.81", dy: 0.02 }],
  caption: "NekIBM Set1sphere: 0.5 mm 유리구가 물에서 정지 상태에서 침강 (Ga 49.1)",
}, [
  ["물리", "Mordant & Pinton의 0.5 mm 유리구, 무차원 d = 1/6, ν 0.005416 (Ga 49.1), ρ_p/ρ_f 2.56, g 9.81, 정지 출발"],
  ["노트북", "1.2 × 3 × 1.2 (7.2d × 18d × 7.2d) 주기, 10 × 25 × 10 요소, N 7 (d/10), 마커 d/17, FWHM 0.12 d, dt 2e-3, 500 스텝 (t = 1)"],
  ["_fine", "같은 격자 N 8, FWHM 0.10 d"],
  ["hires", "1.2 × 10 × 1.2 (60d 기둥), 21 × 180 × 21, N 7 (d/20, NekIBM 격자), 마커 d/20, FWHM 0.10 d, dt 1e-3, 16 MI300A (2.5 s/스텝)"],
  ["기준", "Mordant & Pinton (2000): Re_t 41, u_t 1.332; 표준 항력 ODE (S–N + 부가질량) 1.427 (7 % 높음)"],
], { rowH: 0.7 });
{
  const s = content("V12 결과 — 노트북 −6 %, NekIBM 격자에서 −2.5 % (해상도가 주원인)");
  img(s, "settle_curve.png", 0.5, 1.25, 4.4, 3.57);
  img(s, "settle_slice.png", 5.1, 1.25, 3.9, 3.6);
  grid(s, ["", "v", "Re", "v / u_t"], [
    ["노트북, N 7 (t = 1)", "1.249", "38.4", "0.937"],
    ["_fine, N 8, 0.10 d", "1.244", "38.3", "0.933"],
    ["hires (t = 1 / 1.2)", "1.273 / 1.299", "39.2 / 40.0", "0.956 / 0.975"],
    ["Mordant & Pinton", "1.332", "41", "1"],
  ], 9.2, 1.3, 3.6, [1.5, 0.8, 0.55, 0.75], { fontSize: 9.5, rowH: 0.42 });
  bullets(s, [
    "필터 폭 변화는 0.4 %; NekIBM 격자(d/20)에서 적자가 2.5 %로 줄고 속도는 아직 상승 중",
    "남은 차이는 7.2d 주기 기둥의 방해 효과 크기 (Hasimoto, Stokes 한계 14 %, Re 41에서는 훨씬 작음); 더 넓은 상자로 분리 가능",
    "hires는 pdebug 60분 제한으로 2000 스텝 중 1207에서 종료; 32 랭크는 설정 단계 segfault (16 사용)",
  ], 0.5, 5.2, 12.3, 1.6, 10.5);
  caption(s, "가운데: 구 좌표계 |u − v_b|/|v_b| 단면 (t = 1)", 5.1, 4.9, 3.9);
}

// ---------- V13–V15: Pilobolus 준비 ----------
section = "3. Pilobolus 준비 검증 (V13–V15)";
pres.addSection({ title: section });
{
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: section });
  s.addText("3. Pilobolus 준비 검증 (V13–V15)", { placeholder: "title" });
  s.addText("필터 폭 분리(ppiclF 패치 0006), 처방 운동 물체, 부착 접촉 — 발사 해석에 필요한 세 기능", { placeholder: "body" });
}

// V13
setupSlide("V13 설정 — EL·IBM 필터 폭 분리 (패치 0006, Config::filterWidthIbm)", {
  aspect: 1, faces: { top: "V5 충돌 케이스 (주기 상자)" },
  bodies: [{ cx: 0.55, cy: 0.5, r: 0.15, label: "IBM 마커 폭 σ_IBM (0.03)" }],
  clouds: [{ cx: 0.2, cy: 0.5, w: 0.18, h: 0.18, n: 40, r: 0.03, color: HEX.accent2, label: "EL 투영 폭 σ_EL (0.06)" }],
  caption: "같은 투영 루틴 안에서 입자마다 다른 가우시안 폭 (PPICLF_R_JFILT)",
}, [
  ["변경", "ppiclF 패치 0006: 입자별 가우시안 폭 PPICLF_R_JFILT (실제·고스트 루프), 전역 폭보다 크면 abort"],
  ["인터페이스", "Config::filterWidthIbm; 0이면 EL 폭과 동일"],
  ["케이스 A", "V5 충돌, filter_width_ibm = 0 (동일) — 패치 전·후 비교"],
  ["케이스 B", "V5 충돌, EL 0.06 / IBM 0.03 (IBM 폭 절반)"],
  ["케이스 C", "V11 노트북 Re 100 (마커만) 패치 전·후"],
  ["검사", "구 속도, 운동량 수지, 마커 슬립, 비트 동일성"],
], { rowH: 0.52 });
{
  const s = content("V13 결과 — 동일 폭은 결과 재현, 분리 시 수지 유지·슬립 절반");
  grid(s, ["실행", "구 u (t = 1.2)", "(ΣP − P_0)/P_0", "마커 슬립 (끝)", "비고"], [
    ["동일 폭, 패치 전", "0.27100", "−3.155e-3", "4.843e-3", ""],
    ["동일 폭, 패치 후", "0.27100", "−3.155e-3", "4.844e-3", "다른 자릿수 모두 동일"],
    ["EL 0.06 / IBM 0.03", "0.27176", "−4.542e-3", "6.670e-3", "충돌 중 슬립 1.2e-2 vs 2.0e-2"],
    ["Re 100 C_d (U t/d 10)", "1.2860 전·후", "", "", "비트 동일"],
  ], 0.5, 1.3, 12.3, [2.8, 2.0, 2.2, 2.0, 3.3], { fontSize: 11.5, rowH: 0.45 });
  stat(s, 0.5, 3.7, 3.9, "비트 동일", "동일 폭에서 패치는 inert");
  stat(s, 4.7, 3.7, 3.9, "0.45 %", "분리 시 운동량 수지 (0.32 % → 0.45 %, 같은 지연 기구)");
  stat(s, 8.9, 3.7, 3.9, "2.28 → 1.64", "명시적 물체 갱신의 ρ_b/ρ_f 하한");
  bullets(s, ["Pilobolus 수액 케이스에서 EL 폭 250 µm, IBM 폭 50 µm로 사용 (design.md §7.1)"], 0.5, 5.3, 12.3, 0.8, 12);
}

// V14
setupSlide("V14 설정 — 처방 운동 구 vs 고정 구 (ibm_sphere_moving)", {
  aspect: 0.6, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  bodies: [{ cx: 0.5, cy: 0.78, r: 0.0625, label: "구, y = +5d에서 출발" }],
  arrows: [{ x1: 0.8, y1: 0.68, x2: 0.8, y2: 0.4, label: "v = −u_in (처방)", dy: 0.02 }],
  caption: "V11 상자를 완전 주기로, 유체 정지; 갈릴레이 변환으로 고정 구와 같아야 함",
}, [
  ["기능", "Body::prescribed — 스텝 시작에 UDF가 x(t), v(t)를 주고 동역학은 건너뜀; 유체력은 보고"],
  ["영역·격자", "V11 노트북 상자 (16 × 30 × 16, N 5) 완전 주기, 유체 정지"],
  ["구", "d 0.01, y = +5d에서 −u_in = −1 (y)로 이동, 마커·FWHM은 고정 구와 동일"],
  ["시간", "dt 2e-4, 500 스텝 (U t/d = 10), 4 랭크"],
  ["기준", "같은 격자의 ibm_sphere_re100 (차이: 유입·유출 경계와 시동)"],
  ["측정", "C_d = F_y/(0.5 ρ u_in² π d²/4), U t/d = 8–10 평균 (보고 11회)"],
], { rowH: 0.55 });
{
  const s = content("V14 결과 — 처방 운동 구의 항력이 고정 구와 2.4 % 이내");
  const mv = fs.readFileSync(path.join(__dirname, "moving_cd.txt"), "utf8").trim().split("\n").map((l) => l.split(" ").map(Number));
  const fx = fs.readFileSync(path.join(__dirname, "fixed_cd.txt"), "utf8").trim().split("\n").map((l) => l.split(" ").map(Number));
  s.addChart(pres.ChartType.line, [
    { name: "처방 운동 구 (유체 정지)", labels: mv.map((r) => r[0].toFixed(1)), values: mv.map((r) => Math.abs(r[1])) },
    { name: "고정 구 (유입류)", labels: fx.map((r) => r[0].toFixed(1)), values: fx.map((r) => Math.abs(r[1])) },
  ], {
    x: 0.5, y: 1.3, w: 7.6, h: 4.6, chartColors: [HEX.accent2, HEX.accent1], lineSize: 2, lineDataSymbol: "none",
    showTitle: true, title: "C_d vs U t/d (Re 100, 노트북 격자)", titleFontSize: 12, titleColor: HEX.dk2, titleFontFace: "+mn-lt",
    showLegend: true, legendPos: "b", legendFontSize: 10, legendFontFace: "+mn-lt", legendColor: HEX.dk1,
    catAxisTitle: "U t / d", showCatAxisTitle: true, catAxisTitleFontSize: 10, catAxisTitleColor: HEX.dk1, catAxisLabelFontSize: 9, catAxisLabelColor: HEX.dk1, catAxisLabelFontFace: "+mn-lt", catAxisLabelFrequency: 4,
    valAxisTitle: "C_d", showValAxisTitle: true, valAxisTitleFontSize: 10, valAxisTitleColor: HEX.dk1, valAxisLabelFontSize: 9, valAxisLabelColor: HEX.dk1, valAxisLabelFontFace: "+mn-lt",
    valAxisMinVal: 1.0, valAxisMaxVal: 2.0, valGridLine: { color: "DDE5E8", size: 0.5 }, catGridLine: { style: "none" },
  });
  grid(s, ["U t/d = 8–10", "C_d 평균", "표준편차", "마커 슬립"], [
    ["처방 운동 구", "1.257", "0.029", "3.7–4.0 %"],
    ["고정 구", "1.289", "0.001", "3.9 %"],
  ], 8.4, 1.3, 4.4, [1.5, 1.0, 1.0, 0.9], { fontSize: 11, rowH: 0.4 });
  bullets(s, [
    "±2 %의 스텝별 진동: 움직이는 물체의 마커가 GLL 셀을 가로지르며 명시적 직접 강제가 매 스텝 새 격자 속도를 봄 — 평균되고 해상도로 줄어듦",
    "U t/d = 10의 5 % 하락은 15d 주기 기둥을 돌아 자기 시동 교란에 도달한 것 — 창은 9.8까지",
    "Pilobolus 발사의 처방 v(t) 구동에 이 기능을 사용",
  ], 8.4, 2.8, 4.4, 3.6, 11);
}

// V15
setupSlide("V15 설정 — 부착 접촉: 물체에 포획되는 EL 입자 (el_ibm_adhesion)", {
  aspect: 1, faces: { top: "주기", bottom: "주기", left: "주기", right: "주기" },
  bodies: [{ cx: 0.5, cy: 0.5, r: 0.15, label: "IBM 구 (자유), 붙은 입자가 강체로 동행" }],
  clouds: [{ cx: 0.18, cy: 0.5, w: 0.2, h: 0.2, n: 60, r: 0.03, color: HEX.accent2, label: "EL 200개" }, { cx: 0.36, cy: 0.5, w: 0.03, h: 0.25, n: 14, r: 0.03, color: HEX.accent4 }],
  arrows: [{ x1: 0.3, y1: 0.75, x2: 0.42, y2: 0.75, label: "u = 0.5", dy: 0.02 }],
  caption: "V5와 모든 설정 동일, attach = 1만 켬",
}, [
  ["기능", "Config::attachOnContact — 첫 접촉에서 부착, 물체 속도를 즉시 받음, 잃은 운동량 m_p(u_b − v_p)는 충격량 채널(JIX)로 물체에"],
  ["부착 후", "강체 속도로 이동 (항력·접촉 없음), PlaceMarkers가 강체 위치로 스냅, 부피는 φ_p에 계속 투영, 질량은 물체 관성·무게에 합산"],
  ["EL–EL", "붙은 입자에 작용하는 자유 입자 접촉력은 물체 충격량으로 전달; 대시팟은 강체 속도 사용"],
  ["케이스", "V5 (주기 상자, 구 R 0.15, EL 200개, k_n 500) + attach = 1"],
  ["검사", "포획 수, 반발 수, 최근접 거리, 구 속도, 운동량 수지 (탄성 접촉 V5와 비교)"],
], { rowH: 0.68 });
{
  const s = content("V15 결과 — 포획도 탄성 접촉과 같은 수준으로 운동량을 보존");
  grid(s, ["t = 1.2", "탄성 접촉 (V5)", "부착"], [
    ["붙은 입자", "0", "107 / 200"],
    ["반발 입자", "4", "0"],
    ["최근접 gap/d_p", "0.109", "−0.021 (첫 접촉 겹침)"],
    ["구 속도", "0.271", "0.276"],
    ["(ΣP − P_0)/P_0", "−3.2e-3", "−3.2e-3"],
  ], 0.5, 1.3, 7.0, [2.6, 2.0, 2.4], { fontSize: 12, rowH: 0.45 });
  stat(s, 7.9, 1.3, 2.3, "107 / 200", "포획 (나머지는 빗나감)");
  stat(s, 10.4, 1.3, 2.4, "−3.2e-3", "운동량 수지, 탄성 접촉과 동일");
  stat(s, 7.9, 2.75, 4.9, "2.6x", "붙은 질량 / 구의 초과 질량 — 함께 이동");
  bullets(s, [
    "첫 구현은 18 %를 잃음: 자유 입자가 붙은 입자를 밀 때 힘은 받지만 강체로 움직이는 붙은 입자가 반력을 물체에 넘기지 않았음 → 반력 전달 + 강체 속도 대시팟으로 수정",
    "attach = 0이면 충돌 케이스는 모든 자릿수에서 불변",
    "미시험: 회전하는 물체에 부착, 탈착",
  ], 0.5, 4.3, 12.3, 2.3, 12);
}

// ---------- Pilobolus ----------
section = "4. Pilobolus 포자낭 발사";
pres.addSection({ title: section });
{
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: section });
  s.addText("4. Pilobolus 포자낭 발사", { placeholder: "title" });
  s.addText("Yafetto et al. (2008): d 0.35 mm, 8.6 m/s에 0.22 ms 도달 후 1300 m/s² 감속 — 공기와 수액이 각각 얼마나 감속시키는가", { placeholder: "body" });
}

// P1 setup
setupSlide("P1 설정 — 단상 발사: 처방 v(t)와 자유 비행 (pilobolus_launch)", {
  aspect: 1.2, faces: { top: "원방 u = 0", bottom: "원방 u = 0", left: "원방 u = 0 (노즐 x = 0)", right: "유출 (동압 안정화)" },
  band: { x: 0, y: 0.44, w: 1, h: 0.12, label: "세밀 띠 |y|,|z| < 0.6 mm (0.2 mm 요소, N 5: 40 µm ≈ d/9)" },
  bodies: [{ cx: 0.19, cy: 0.5, r: 0.03 }],
  arrows: [{ x1: 0.27, y1: 0.5, x2: 0.6, y2: 0.5, label: "포자낭 d 0.35 mm, v(t): 0.22 ms에 8.6 m/s, 이후 −1300 m/s²", dy: 0.62, lw: 4.2 }],
  caption: "x −2~10 mm, y·z ±5 mm, 11 760 요소 (hires: 105 800, 0.08 mm·N 7 = 11 µm); 8 랭크 (macOS 2 GB/랭크 한계)",
}, [
  ["공기", "ρ 1.2, μ 1.8e-5; dt 1 µs (hires 0.3 µs), 1.2 ms"],
  ["포자낭", "등가 구 d 350 µm, ρ_b 1150; 마커 25 µm, FWHM 50 µm (hires 8 / 16 µm)"],
  ["처방 (mode 0)", "x(t), v(t)를 전 구간 강제; 보고: 마커 힘 합 F_m, 물체 힘 F = F_m + ρ_f V a, C_d, a_f = F/(ρ_b V)"],
  ["자유 비행 (mode 1)", "t_rise = 0.22 ms에 8.6 m/s로 놓음 → 유체력·중력만으로 감속"],
  ["측정 (Yafetto 2008)", "v_max 8.6 m/s, t_rise 0.22 ms, 활공 감속 ~1300 m/s², Re ≈ 200"],
  ["실행", "노트북 8 랭크 CPU, 3.4–6.8 s/스텝 (각 1200 스텝); hires는 Tuolumne 입력 준비됨"],
], { rowH: 0.6 });
{
  const s = content("P1 결과 — 공기 저항은 측정 감속의 9 %뿐: 포자낭은 거의 감속 없이 활공");
  img(s, "pil_top.png", 0.5, 1.25, 6.6, 4.5);
  grid(s, ["활공 구간 (t > 0.22 ms)", "처방 v(t)", "자유 비행", "측정"], [
    ["1.2 ms 속도", "7.33 m/s (강제)", "8.48 m/s", "~7.3 m/s"],
    ["유체력", "2.4–3.3 µN", "약간 큼 (속도 유지)", ""],
    ["감속", "a_f 111 m/s² (93–141)", "125 m/s²", "1300 m/s²"],
    ["C_d", "0.785", "0.70–0.75", "S–N 0.832 (Re 185)"],
  ], 7.3, 1.3, 5.5, [1.9, 1.3, 1.2, 1.1], { fontSize: 10, rowH: 0.4 });
  stat(s, 7.3, 3.5, 2.65, "9 %", "공기가 설명하는 측정 감속 비율");
  stat(s, 10.15, 3.5, 2.65, "1190 m/s²", "수액 쪽에서 나와야 하는 감속");
  bullets(s, [
    "두 모드가 같은 답: 건조 포자낭은 0.96 ms 동안 0.12 m/s만 잃음",
    "가속 구간 마커 힘 −5.3 µN 중 2.1 µN은 껍질 안 유체의 관성 (ρ_f V a); 물체 힘 최대 −4.2 µN",
    "힘의 ±0.2 µN 톱니(주기 27 µs)는 0.2 mm 요소를 7.5 m/s로 건너는 격자 교차 잡음 (V14와 같은 기구)",
  ], 7.3, 5.0, 5.5, 1.9, 10.5);
  caption(s, "scripts/pilobolus_launch.py → docs/img/pilobolus_launch.png (위 두 줄)", 0.5, 5.85, 6.6);
}

// P2 setup
setupSlide("P2 설정 — 수액 액적 주입과 부착 (pilobolus_launch/sap)", {
  aspect: 1.2, faces: { top: "원방 u = 0", bottom: "원방 u = 0", left: "노즐 r 75 µm (x = 0)", right: "유출" },
  band: { x: 0, y: 0.44, w: 1, h: 0.12, label: "세밀 띠 (P1과 동일 격자)" },
  bodies: [{ cx: 0.3, cy: 0.5, r: 0.03 }],
  clouds: [{ cx: 0.14, cy: 0.5, w: 0.2, h: 0.05, n: 30, r: 0.045, color: HEX.accent2 }, { cx: 0.14, cy: 0.5, w: 0.22, h: 0.08, n: 60, r: 0.02, color: HEX.accent5 }],
  arrows: [{ x1: 0.02, y1: 0.36, x2: 0.12, y2: 0.36, label: "수액 제트 12 m/s, 0.3 ms", dy: 0.72, lw: 2.4 }],
  caption: "P1의 단상 케이스 + EL 수액 액적 (주액적 120 µm, 위성 18 µm): 노즐 원반에서 매 스텝 주입, 포자낭(처방 v(t))에 접촉 시 부착 (V15)",
}, [
  ["수액 (Yafetto 2008)", "총 30 nL (10–50 nL), ρ 1005; 주액적 d 120 µm (부피 90 %, 30개), 위성 18 µm (10 %, ~980개)"],
  ["주입", "노즐 반지름 75 µm, v_jet 12 m/s, 횡 분산 0.5 m/s, t < 0.3 ms 동안 sap_volume·dt/t_inj 씩 (소수 개수 이월)"],
  ["필터 폭", "EL 투영 250 µm (액적 간격 몇 개), IBM 마커 50 µm — V13의 분리"],
  ["접촉·부착", "attach = 1 (첫 접촉에서 부착, V15), k_n 2000, 서브스텝 5"],
  ["보고 추가", "영역 내 액적 수, 붙은 수·질량, 수액 운동량; 스냅샷 5 스텝, 체크포인트 25 스텝 (영상)"],
  ["실행", "노트북 8 랭크, 3–23 s/스텝 (주입 중 느림), 1200 스텝 ≈ 2.5 h"],
], { rowH: 0.6 });
{
  const s = content("P2 결과 — 수액은 질량의 90 %만큼 붙지만 감속에는 기여하지 못함");
  img(s, "pil_sap.png", 0.5, 1.25, 6.4, 2.18);
  img(s, "sap_f060.png", 0.5, 3.55, 3.1, 1.16);
  img(s, "sap_f200.png", 3.8, 3.55, 3.1, 1.16);
  grid(s, ["항목", "값"], [
    ["주입 액적", "1011개 (주 30, 위성 ~980), 0.3 ms"],
    ["붙은 액적", "775개, 0.86 ms에 포화; 끝에 1005개 잔류"],
    ["붙은 질량", "23.2 ng = 수액의 77 %, 포자낭 질량의 90 %"],
    ["활공 유체력", "1.8 µN vs 건조 2.9 µN (−36 %)"],
    ["C_d", "0.4–0.6 vs 건조 0.79"],
    ["공기만의 감속 (붙은 질량 포함)", "~40 m/s² (측정 1300의 3 %)"],
  ], 7.2, 1.3, 5.6, [2.4, 3.2], { fontSize: 10, rowH: 0.36 });
  bullets(s, [
    "수액 30 nL(30 ng)이 포자낭(26 ng)보다 무거워 붙은 질량이 관성을 키움 → 감속 감소; 12 m/s 액적의 충돌 운동량은 앞 방향",
    "뒤따르는 제트가 포자낭 뒤 공기를 끌고 가 상대속도가 줄고 항력이 36 % 감소",
    "결론: 붙은 수액의 질량·항력으로는 측정 감속을 설명 못 함 → 수액 꼬리의 장력(테더) 또는 수액 원반의 유효 단면이 다음 모델 (미구현)",
  ], 7.2, 4.6, 5.6, 2.3, 10);
  caption(s, "위: 액적 수·붙은 질량·수액 운동량. 아래: ParaView 영상 프레임 t = 0.30, 1.00 ms (회색 구 포자낭, 주황 액적, φ_p 단면)", 0.5, 5.0, 6.5);
}

// ---------- summary ----------
section = "5. 요약";
pres.addSection({ title: section });
{
  const s = content("요약 — 검증 상태와 남은 일");
  grid(s, ["영역", "결론"], [
    ["EL 결합 (V1–V9)", "배관 반올림 수준; 보존 수지 1e-3; 두 독립 경로 1 % 이내; Nek5000과 1.2 % 이내; Bickley 3.5 %; TAB 개시 곡선 위"],
    ["IBM (V10–V12)", "Re 10 기준값을 폭으로 감쌈; Re 100 C_d +18 → +10 % (해상도); 침강 −6 → −2.5 % (NekIBM 격자)"],
    ["Pilobolus 준비 (V13–V15)", "필터 폭 분리 inert·수지 유지; 처방 운동 구 2.4 % 이내; 부착 접촉 수지 동일"],
    ["Pilobolus 발사 (P1, P2)", "공기 9 %, 붙은 수액 3 % → 측정 감속 1300 m/s²의 나머지는 수액 꼬리 장력 등 다른 기구"],
  ], 0.5, 1.3, 12.3, [2.8, 9.5], { fontSize: 11.5, rowH: 0.62 });
  bullets(s, [
    "Tuolumne 대기: Re 300 고정 구 (Johnson & Patel C_d 0.656, St 0.137), Pilobolus hires 발사 (105 800 요소, N 7)",
    "모델: 수액 꼬리 장력(테더) 구현 여부; 주액적(120 µm, 격자 3칸)은 검증된 점입자 영역 밖 — 40 µm 액적 민감도 실행 제안",
    "미해결: 32 랭크 설정 segfault, TAB We ≥ 200 재분열 불안정, 디바이스 경로의 GPU 백엔드 실행, NekIBM 자체 실행 (Linux/Intel 필요)",
  ], 0.5, 5.1, 12.3, 1.9, 12);
}

(async () => {
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log("wrote", OUT);
})();
