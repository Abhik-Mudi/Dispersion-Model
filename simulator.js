// ============================================================================
// 1. DISPERSION COEFFICIENTS (PASQUILL-GIFFORD CORRELATIONS)
// ============================================================================

function getPlumeSigmas(x_m, stabilityClass = "D", terrain = "rural") {
  const x = Math.max(Number(x_m), 0.1);
  const sc = String(stabilityClass).toUpperCase();
  const terr = String(terrain).toLowerCase();

  let sy = 0.0;
  let sz = 0.0;

  if (terr === "rural") {
    if (sc === "A") {
      sy = 0.22 * x * Math.pow(1.0 + 0.0001 * x, -0.5);
      sz = 0.20 * x;
    } else if (sc === "B") {
      sy = 0.16 * x * Math.pow(1.0 + 0.0001 * x, -0.5);
      sz = 0.12 * x;
    } else if (sc === "C") {
      sy = 0.11 * x * Math.pow(1.0 + 0.0001 * x, -0.5);
      sz = 0.08 * x * Math.pow(1.0 + 0.0002 * x, -0.5);
    } else if (sc === "D") {
      sy = 0.08 * x * Math.pow(1.0 + 0.0001 * x, -0.5);
      sz = 0.06 * x * Math.pow(1.0 + 0.0015 * x, -0.5);
    } else if (sc === "E") {
      sy = 0.06 * x * Math.pow(1.0 + 0.0001 * x, -0.5);
      sz = 0.03 * x * Math.pow(1.0 + 0.0003 * x, -1.0);
    } else if (sc === "F") {
      sy = 0.04 * x * Math.pow(1.0 + 0.0001 * x, -0.5);
      sz = 0.016 * x * Math.pow(1.0 + 0.0003 * x, -1.0);
    } else {
      throw new Error(`Invalid stability class '${sc}'. Choose A to F.`);
    }
  } else if (terr === "urban") {
    if (sc === "A" || sc === "B") {
      sy = 0.32 * x * Math.pow(1.0 + 0.0004 * x, -0.5);
      sz = 0.24 * x * Math.pow(1.0 + 0.001 * x, 0.5);
    } else if (sc === "C") {
      sy = 0.22 * x * Math.pow(1.0 + 0.0004 * x, -0.5);
      sz = 0.20 * x;
    } else if (sc === "D") {
      sy = 0.16 * x * Math.pow(1.0 + 0.0004 * x, -0.5);
      sz = 0.14 * x * Math.pow(1.0 + 0.0003 * x, -0.5);
    } else if (sc === "E" || sc === "F") {
      sy = 0.11 * x * Math.pow(1.0 + 0.0004 * x, -0.5);
      sz = 0.08 * x * Math.pow(1.0 + 0.0015 * x, -0.5);
    } else {
      throw new Error(`Invalid stability class '${sc}'. Choose A to F.`);
    }
  } else {
    throw new Error("Terrain must be 'rural' or 'urban'.");
  }

  return { sy, sz };
}

const PUFF_COEFFS = {
  A: [0.18, 0.92, 0.60, 0.75],
  B: [0.14, 0.92, 0.53, 0.73],
  C: [0.10, 0.92, 0.34, 0.71],
  D: [0.06, 0.92, 0.15, 0.70],
  E: [0.04, 0.92, 0.10, 0.65],
  F: [0.02, 0.89, 0.05, 0.61],
};

function getPuffSigmas(x_m, stabilityClass = "D") {
  const x = Math.max(Number(x_m), 0.1);
  const sc = String(stabilityClass).toUpperCase();
  const coeffs = PUFF_COEFFS[sc] || PUFF_COEFFS.D;
  const [ay, by, cz, dz] = coeffs;

  const sy = ay * Math.pow(x, by);
  const sz = cz * Math.pow(x, dz);
  const sx = sy; // Isotropic horizontal assumption
  return { sx, sy, sz };
}

// ============================================================================
// 2. CONCENTRATION EVALUATION C(x, y, z, t)
// ============================================================================

function computeConcentration(
  modelType,
  sourceStrength,
  u_m_s,
  Hr_m,
  x_m,
  y_m,
  z_m,
  t_sec,
  stabilityClass = "D",
  terrain = "rural"
) {
  const X = Number(x_m);
  const Y = Number(y_m);
  const Z = Number(z_m);
  const u = Math.max(Number(u_m_s), 0.001);
  const Hr = Number(Hr_m);
  const t = Math.max(Number(t_sec), 0.0);
  const Q = Number(sourceStrength);
  const type = String(modelType).toLowerCase();

  if (type === "puff") {
    const travelDist = Math.max(u * t, 0.1);
    const evalDist = Math.max(X, travelDist);
    const { sx, sy, sz } = getPuffSigmas(evalDist, stabilityClass);

    const denom = Math.pow(2.0 * Math.PI, 1.5) * sx * sy * sz;
    const exp_x = Math.exp(-0.5 * Math.pow((X - u * t) / sx, 2));
    const exp_y = Math.exp(-0.5 * Math.pow(Y / sy, 2));
    const exp_z =
      Math.exp(-0.5 * Math.pow((Z - Hr) / sz, 2)) +
      Math.exp(-0.5 * Math.pow((Z + Hr) / sz, 2));

    return (Q / denom) * exp_x * exp_y * exp_z;
  } else if (type === "plume") {
    const { sy, sz } = getPlumeSigmas(X, stabilityClass, terrain);
    const denom = 2.0 * Math.PI * sy * sz * u;
    const exp_y = Math.exp(-0.5 * Math.pow(Y / sy, 2));
    const exp_z =
      Math.exp(-0.5 * Math.pow((Z - Hr) / sz, 2)) +
      Math.exp(-0.5 * Math.pow((Z + Hr) / sz, 2));
    const C_steady = (Q / denom) * exp_y * exp_z;

    const plumeFrontX = u * t;
    return X > 0.0 && X <= plumeFrontX ? C_steady : 0.0;
  } else {
    throw new Error("modelType must be 'plume' or 'puff'.");
  }
}

// ============================================================================
// 3. CHEMICAL PRESETS
// ============================================================================

const CHEMICAL_PRESETS = {
  Chlorine: { mw: 70.9, defaultThreshold: 3.0, desc: "Chlorine gas (Cl₂)" },
  "Sulfur Dioxide": { mw: 64.07, defaultThreshold: 0.5, desc: "Sulfur dioxide (SO₂)" },
  Ammonia: { mw: 17.03, defaultThreshold: 150.0, desc: "Anhydrous ammonia (NH₃)" },
  "Hydrogen Sulfide": { mw: 34.08, defaultThreshold: 30.0, desc: "Hydrogen sulfide (H₂S)" },
  Benzene: { mw: 78.11, defaultThreshold: 150.0, desc: "Benzene vapor (C₆H₆)" },
  "Carbon Monoxide": { mw: 28.01, defaultThreshold: 350.0, desc: "Carbon monoxide (CO)" },
  Methane: { mw: 16.04, defaultThreshold: 50000.0, desc: "Methane gas (CH₄)" },
  Custom: { mw: 70.9, defaultThreshold: 3.0, desc: "User-defined chemical" },
};

// ============================================================================
// 3.5. PASQUILL-GIFFORD TABLE 5-1 STABILITY PREDICTION ENGINE
// ============================================================================

const TABLE_5_1 = {
  rows: [
    { bracket: "< 2", min: 0.0, max: 2.0, name: "Calm / Light air (< 2 m/s)" },
    { bracket: "2 – 3", min: 2.0, max: 3.0, name: "Light breeze (2–3 m/s)" },
    { bracket: "3 – 4", min: 3.0, max: 4.0, name: "Gentle breeze (3–4 m/s)" },
    { bracket: "4 – 6", min: 4.0, max: 6.0, name: "Moderate breeze (4–6 m/s)" },
    { bracket: "> 6", min: 6.0, max: Infinity, name: "Strong breeze (> 6 m/s)" },
  ],
  grid: [
    // Row 0: < 2 m/s
    { day_strong: "A", day_moderate: "A-B", day_slight: "B", night_overcast: "F", night_clear: "F" },
    // Row 1: 2 - 3 m/s
    { day_strong: "A-B", day_moderate: "B", day_slight: "C", night_overcast: "E", night_clear: "F" },
    // Row 2: 3 - 4 m/s
    { day_strong: "B", day_moderate: "B-C", day_slight: "C", night_overcast: "D", night_clear: "E" },
    // Row 3: 4 - 6 m/s
    { day_strong: "C", day_moderate: "C-D", day_slight: "D", night_overcast: "D", night_clear: "D" },
    // Row 4: > 6 m/s
    { day_strong: "C", day_moderate: "D", day_slight: "D", night_overcast: "D", night_clear: "D" },
  ],
};

const STABILITY_CLASS_NAMES = {
  A: { name: "Extremely Unstable", desc: "Strong convective mixing, light wind & intense solar radiation" },
  B: { name: "Moderately Unstable", desc: "Moderate convective mixing, warm sunny surface conditions" },
  C: { name: "Slightly Unstable", desc: "Slight convective mixing, partly cloudy daylight or low sun angle" },
  D: { name: "Neutral Stability", desc: "Neutral conditions; overcast day/night or higher wind speeds" },
  E: { name: "Slightly Stable", desc: "Weak ground inversion, light wind and nighttime conditions" },
  F: { name: "Moderately Stable", desc: "Strong ground temperature inversion, calm clear night sky" },
};

const ENV_CONDITION_NAMES = {
  day_strong: "Daytime Strong Insolation (Sunny midday)",
  day_moderate: "Daytime Moderate Insolation (Moderate sun)",
  day_slight: "Daytime Slight Insolation (Low sun / overcast)",
  night_overcast: "Nighttime Thin Overcast or > 4/8 Low Cloud",
  night_clear: "Nighttime Clear Sky (≤ 3/8 Cloudiness)",
};

function getWindRowIndex(windSpeed) {
  const u = Math.max(Number(windSpeed) || 0.0, 0.0);
  if (u < 2.0) return 0;
  if (u < 3.0) return 1;
  if (u < 4.0) return 2;
  if (u <= 6.0) return 3;
  return 4;
}

function evaluateStabilityPrediction(windSpeed, envCond, preferredSplit = null) {
  const rowIndex = getWindRowIndex(windSpeed);
  const rowData = TABLE_5_1.grid[rowIndex];
  const rawClass = rowData[envCond] || "D";

  let appliedClass = rawClass;
  let isSplit = false;
  let splitOptions = [];

  if (rawClass.includes("-")) {
    isSplit = true;
    splitOptions = rawClass.split("-");
    if (preferredSplit && splitOptions.includes(preferredSplit)) {
      appliedClass = preferredSplit;
    } else {
      // Default to the second class for conservative safety (e.g. B in A-B, C in B-C, D in C-D)
      appliedClass = splitOptions[1];
    }
  }

  return {
    rowIndex,
    colKey: envCond,
    rawClass,
    appliedClass,
    isSplit,
    splitOptions,
    rowMeta: TABLE_5_1.rows[rowIndex],
    classMeta: STABILITY_CLASS_NAMES[appliedClass] || STABILITY_CLASS_NAMES.D,
    condName: ENV_CONDITION_NAMES[envCond] || envCond,
  };
}

function setStabilityMode(mode) {
  AppState.stabilityMode = mode;
  const toggleCheckbox = document.getElementById("stabilityModeToggle");
  const btnDirect = document.getElementById("btnModeDirect");
  const btnAuto = document.getElementById("btnModeAuto");
  const directSec = document.getElementById("directStabilitySection");
  const autoSec = document.getElementById("autoStabilitySection");

  if (mode === "auto") {
    if (toggleCheckbox) toggleCheckbox.checked = true;
    if (btnDirect) btnDirect.classList.remove("active");
    if (btnAuto) btnAuto.classList.add("active");
    if (directSec) directSec.style.display = "none";
    if (autoSec) autoSec.style.display = "block";
    updateEnvironmentalPrediction();
  } else {
    if (toggleCheckbox) toggleCheckbox.checked = false;
    if (btnDirect) btnDirect.classList.add("active");
    if (btnAuto) btnAuto.classList.remove("active");
    if (directSec) directSec.style.display = "block";
    if (autoSec) autoSec.style.display = "none";
    updateDirectStabilityHint();
  }
}

function updateEnvironmentalPrediction() {
  const uWindInput = document.getElementById("uWind");
  const windVal = parseFloat(uWindInput?.value) || 3.0;
  const selectedCondRadio = document.querySelector('input[name="envCondition"]:checked');
  const envCond = selectedCondRadio?.value || AppState.envCondition || "night_overcast";
  AppState.envCondition = envCond;

  const pred = evaluateStabilityPrediction(windVal, envCond, AppState.preferredSplit);

  // Update Wind Speed Display and Stepper
  const envWindDisplay = document.getElementById("envWindDisplay");
  if (envWindDisplay) envWindDisplay.innerText = `${windVal.toFixed(1)} m/s`;

  const envWindBracketPill = document.getElementById("envWindBracketPill");
  if (envWindBracketPill) envWindBracketPill.innerText = pred.rowMeta.name;

  document.querySelectorAll("#envWindStepper .wind-step").forEach((step, idx) => {
    step.classList.toggle("active", idx === pred.rowIndex);
  });

  // Update Prediction Pill and Titles
  const predClassPill = document.getElementById("predClassPill");
  if (predClassPill) {
    predClassPill.innerText = `Class ${pred.rawClass}`;
    predClassPill.className = `pred-pill ${pred.appliedClass}`;
  }

  const predClassTitle = document.getElementById("predClassTitle");
  if (predClassTitle) {
    predClassTitle.innerText = `${pred.classMeta.name} (Applying Class ${pred.appliedClass})`;
  }

  const predClassSub = document.getElementById("predClassSub");
  if (predClassSub) {
    predClassSub.innerText = pred.classMeta.desc;
  }

  // Update Formula Explanation
  const predFormula = document.getElementById("predFormulaExplanation");
  if (predFormula) {
    predFormula.innerHTML = `
      Match: Wind <strong>${windVal.toFixed(1)} m/s</strong> (${pred.rowMeta.bracket} m/s bracket) &times; <strong>${pred.condName}</strong> &rarr; Predicts <strong>Class ${pred.rawClass}</strong> (Using <strong>Class ${pred.appliedClass}</strong>).
    `;
  }

  // Split Class Resolution (e.g. A-B, B-C, C-D)
  const splitSelector = document.getElementById("predSplitSelector");
  const splitBtnsWrap = document.getElementById("predSplitBtns");
  const splitText = document.getElementById("predSplitText");
  if (splitSelector && splitBtnsWrap) {
    if (pred.isSplit) {
      splitSelector.style.display = "flex";
      if (splitText) splitText.innerText = pred.rawClass;
      splitBtnsWrap.innerHTML = pred.splitOptions
        .map(
          (opt) =>
            `<button type="button" class="split-btn ${opt === pred.appliedClass ? "active" : ""}" data-split="${opt}">
              Use Class ${opt} ${opt === pred.splitOptions[1] ? "(Standard)" : ""}
            </button>`
        )
        .join("");

      splitBtnsWrap.querySelectorAll(".split-btn").forEach((btn) => {
        btn.addEventListener("click", () => {
          AppState.preferredSplit = btn.getAttribute("data-split");
          updateEnvironmentalPrediction();
          if (AppState.autoCalculate) runSimulation();
        });
      });
    } else {
      splitSelector.style.display = "none";
      AppState.preferredSplit = null;
    }
  }

  // Auto-check the stability radio button
  const matchingRadio = document.querySelector(`input[name="stability"][value="${pred.appliedClass}"]`);
  if (matchingRadio) {
    matchingRadio._autoSetting = true;
    matchingRadio.checked = true;
    setTimeout(() => {
      matchingRadio._autoSetting = false;
    }, 50);
  }

  // Highlight active row, column, and cell in Table 5-1
  updateMatrixTableHighlight(pred.rowIndex, pred.colKey);

  // Update hint in direct section too
  updateDirectStabilityHint(pred.appliedClass);
}

function updateMatrixTableHighlight(rowIndex, colKey) {
  const table = document.getElementById("pgMatrixTable");
  if (!table) return;

  // Clear previous highlights
  table.querySelectorAll(".highlight-row").forEach((el) => el.classList.remove("highlight-row"));
  table.querySelectorAll(".highlight-col").forEach((el) => el.classList.remove("highlight-col"));
  table.querySelectorAll(".highlight-active").forEach((el) => el.classList.remove("highlight-active"));

  // Highlight active column header
  const th = table.querySelector(`th[data-col="${colKey}"]`);
  if (th) th.classList.add("highlight-col");

  // Highlight active row
  const tr = table.querySelector(`tr[data-row="${rowIndex}"]`);
  if (tr) {
    tr.classList.add("highlight-row");
    // Highlight intersecting cell
    const cell = tr.querySelector(`td[data-col="${colKey}"]`);
    if (cell) cell.classList.add("highlight-active");
  }
}

function updateDirectStabilityHint(customClass = null) {
  const hintEl = document.getElementById("directStabilityHint");
  if (!hintEl) return;
  const activeClass = customClass || document.querySelector('input[name="stability"]:checked')?.value || "D";
  const meta = STABILITY_CLASS_NAMES[activeClass] || STABILITY_CLASS_NAMES.D;
  hintEl.innerHTML = `Selected: <strong>Class ${activeClass}</strong> — ${meta.name} (${meta.desc})`;
}

// ============================================================================
// 4. APPLICATION STATE & CONTROLLER
// ============================================================================

const AppState = {
  activeEngine: "client",
  backendAvailable: false,
  autoCalculate: true,
  lastResult: null,
  activeTab: "contour", // "contour" | "breakthrough" | "profiles"
  profileSubTab: "crosswind", // "crosswind" | "downwind"
  hoverCoord: null,
  stabilityMode: "direct", // "direct" | "auto"
  envCondition: "night_overcast",
  preferredSplit: null,
};

// Helper: read form data
function getFormData() {
  const chemSelect = document.getElementById("chemSelect");
  const chemNameInput = document.getElementById("chemName");
  const mwInput = document.getElementById("mwInput");
  const modelType = document.querySelector('input[name="modelType"]:checked')?.value || "puff";
  const sourceStrengthInput = document.getElementById("sourceStrength");
  const sourceUnit = document.getElementById("sourceUnit")?.value || "kg";
  const uWindInput = document.getElementById("uWind");
  const stabilityInput = document.querySelector('input[name="stability"]:checked')?.value || "D";
  const terrainInput = document.querySelector('input[name="terrain"]:checked')?.value || "rural";
  const hrInput = document.getElementById("hrInput");
  const thresholdInput = document.getElementById("thresholdInput");
  const xQueryInput = document.getElementById("xQuery");
  const yQueryInput = document.getElementById("yQuery");
  const zQueryInput = document.getElementById("zQuery");
  const tQueryInput = document.getElementById("tQuery");

  let rawStrength = parseFloat(sourceStrengthInput.value) || 0.0;
  // Convert unit to base kg or kg/s
  let strength = rawStrength;
  if (sourceUnit === "g" || sourceUnit === "g/s") {
    strength = rawStrength / 1000.0;
  }

  return {
    chem_name: chemNameInput.value.trim() || chemSelect.value,
    mw_g_mol: Math.max(parseFloat(mwInput.value) || 70.9, 0.1),
    model_type: modelType,
    source_strength: Math.max(strength, 0.0),
    raw_strength: rawStrength,
    source_unit: sourceUnit,
    u_wind: Math.max(parseFloat(uWindInput.value) || 3.0, 0.01),
    stability: stabilityInput,
    terrain: terrainInput,
    hr_m: Math.max(parseFloat(hrInput.value) || 0.0, 0.0),
    c_target_ppm: Math.max(parseFloat(thresholdInput.value) || 3.0, 1e-6),
    x_query: Math.max(parseFloat(xQueryInput.value) || 300.0, 0.0),
    y_query: parseFloat(yQueryInput.value) || 0.0,
    z_query: Math.max(parseFloat(zQueryInput.value) || 0.0, 0.0),
    t_query: Math.max(parseFloat(tQueryInput.value) || 100.0, 0.0),
  };
}

// Run simulation calculation
async function runSimulation() {
  const params = getFormData();
  const runBtn = document.getElementById("runBtn");
  if (runBtn) {
    runBtn.classList.add("running");
    runBtn.innerText = "Simulating...";
  }

  try {
    let result = null;

    if (!result || !result.success) {
      // Run client calculation
      result = simulateClient(params);
    }

    AppState.lastResult = result;
    renderResults(result);
  } catch (error) {
    console.error("Simulation error:", error);
    alert("Simulation calculation error: " + error.message);
  } finally {
    if (runBtn) {
      runBtn.classList.remove("running");
      runBtn.innerHTML = `<span>⚡</span> Run Simulation`;
    }
  }
}

// Client-side simulation pipeline
function simulateClient(params) {
  const {
    chem_name,
    mw_g_mol,
    model_type,
    source_strength,
    u_wind,
    stability,
    terrain,
    hr_m,
    c_target_ppm,
    x_query,
    y_query,
    z_query,
    t_query,
  } = params;

  // 1. Point concentration
  const c_point_kg_m3 = computeConcentration(
    model_type,
    source_strength,
    u_wind,
    hr_m,
    x_query,
    y_query,
    z_query,
    t_query,
    stability,
    terrain
  );
  const c_point_mg_m3 = c_point_kg_m3 * 1e6;
  const c_point_ppm = (c_point_kg_m3 * 24.45 * 1e6) / mw_g_mol;
  const exceeds = c_point_ppm >= c_target_ppm;
  const hazardRatio = c_target_ppm > 0 ? (c_point_ppm / c_target_ppm) * 100 : 0;

  const t_arrival = u_wind > 0 ? x_query / u_wind : 0.0;
  const advection_front_x = u_wind * t_query;
  const has_arrived = advection_front_x >= x_query;

  let sigmas;
  if (model_type === "plume") {
    const { sy, sz } = getPlumeSigmas(x_query, stability, terrain);
    sigmas = { sx: sy, sy, sz };
  } else {
    const travelDist = Math.max(u_wind * t_query, 0.1);
    const evalDist = Math.max(x_query, travelDist);
    sigmas = getPuffSigmas(evalDist, stability);
  }

  // 2. Temporal history C(t) at receptor point
  const t_max = Math.max(t_query * 2.0, t_arrival * 2.0, 60.0);
  const t_steps = 150;
  const t_vec = [];
  const c_time_ppm = [];
  let peak_ppm = 0.0;
  let t_peak = 0.0;

  for (let i = 0; i <= t_steps; i++) {
    const t_val = 0.1 + (i / t_steps) * (t_max - 0.1);
    t_vec.push(t_val);
    const c_kg = computeConcentration(
      model_type,
      source_strength,
      u_wind,
      hr_m,
      x_query,
      y_query,
      z_query,
      t_val,
      stability,
      terrain
    );
    const ppm = (c_kg * 24.45 * 1e6) / mw_g_mol;
    c_time_ppm.push(ppm);
    if (ppm > peak_ppm) {
      peak_ppm = ppm;
      t_peak = t_val;
    }
  }

  // 3. Crosswind Profile C(y) at x_query, z_query, t_query
  const y_max = Math.max(Math.abs(y_query) * 2.5, 60.0, sigmas.sy * 3.5);
  const y_steps = 80;
  const y_vec = [];
  const c_y_ppm = [];
  for (let i = 0; i <= y_steps; i++) {
    const y_val = -y_max + (i / y_steps) * (2 * y_max);
    y_vec.push(y_val);
    const c_kg = computeConcentration(
      model_type,
      source_strength,
      u_wind,
      hr_m,
      x_query,
      y_val,
      z_query,
      t_query,
      stability,
      terrain
    );
    c_y_ppm.push((c_kg * 24.45 * 1e6) / mw_g_mol);
  }

  // 4. Downwind Profile C(x) at y=0, z_query, t_query
  const x_max_domain = Math.max(x_query * 1.8, advection_front_x * 1.6, 350.0);
  const x_steps = 90;
  const x_vec = [];
  const c_x_ppm = [];
  for (let i = 0; i <= x_steps; i++) {
    const x_val = 1.0 + (i / x_steps) * (x_max_domain - 1.0);
    x_vec.push(x_val);
    const c_kg = computeConcentration(
      model_type,
      source_strength,
      u_wind,
      hr_m,
      x_val,
      0.0,
      z_query,
      t_query,
      stability,
      terrain
    );
    c_x_ppm.push((c_kg * 24.45 * 1e6) / mw_g_mol);
  }

  // 5. 2D Spatial Grid C(x, y) at z_query, t_query
  const grid_nx = 60;
  const grid_ny = 60;
  const grid_x = [];
  const grid_y = [];
  const c_matrix = [];
  let max_ppm = 0.0;

  for (let ix = 0; ix <= grid_nx; ix++) {
    grid_x.push(1.0 + (ix / grid_nx) * (x_max_domain - 1.0));
  }
  for (let iy = 0; iy <= grid_ny; iy++) {
    grid_y.push(-y_max + (iy / grid_ny) * (2 * y_max));
  }

  for (let iy = 0; iy <= grid_ny; iy++) {
    const row = [];
    const y_val = grid_y[iy];
    for (let ix = 0; ix <= grid_nx; ix++) {
      const x_val = grid_x[ix];
      const c_kg = computeConcentration(
        model_type,
        source_strength,
        u_wind,
        hr_m,
        x_val,
        y_val,
        z_query,
        t_query,
        stability,
        terrain
      );
      const ppm = (c_kg * 24.45 * 1e6) / mw_g_mol;
      row.push(ppm);
      if (ppm > max_ppm) max_ppm = ppm;
    }
    c_matrix.push(row);
  }

  return {
    success: true,
    inputs: params,
    point: {
      c_kg_m3: c_point_kg_m3,
      c_mg_m3: c_point_mg_m3,
      c_ppm: c_point_ppm,
      exceeds,
      hazard_ratio_percent: hazardRatio,
      t_arrival_sec: t_arrival,
      advection_front_m: advection_front_x,
      has_arrived,
      sigmas,
      peak_ppm,
      t_peak_sec: t_peak,
    },
    temporal: {
      t_vec,
      c_vec: c_time_ppm,
    },
    crosswind: {
      y_vec,
      c_vec: c_y_ppm,
    },
    downwind: {
      x_vec,
      c_vec: c_x_ppm,
    },
    grid: {
      x_vec: grid_x,
      y_vec: grid_y,
      c_matrix,
      max_ppm,
    },
  };
}

// ============================================================================
// 5. UI RENDERING & VISUALIZATIONS
// ============================================================================

function formatSci(num, decimals = 4) {
  if (num === 0) return "0.0000";
  if (Math.abs(num) >= 0.01 && Math.abs(num) < 10000) {
    return num.toFixed(decimals);
  }
  return num.toExponential(decimals);
}

function renderResults(res) {
  if (!res || !res.point) return;
  const p = res.point;
  const inp = res.inputs;

  // 1. Hazard Alert Banner
  const alertBanner = document.getElementById("hazardAlert");
  if (alertBanner) {
    if (p.exceeds) {
      alertBanner.className = "alert-banner danger";
      alertBanner.innerHTML = `
        <div class="alert-icon">⚠️</div>
        <div class="alert-content">
          <h4>DANGER: EXCEEDS HAZARD THRESHOLD</h4>
          <p>Concentration at monitoring station is <strong>${formatSci(p.c_ppm, 4)} ppm</strong>, exceeding the safety limit of <strong>${inp.c_target_ppm.toFixed(2)} ppm</strong> (${p.hazard_ratio_percent.toFixed(1)}% of threshold).</p>
        </div>
      `;
    } else if (p.c_ppm > 0) {
      alertBanner.className = "alert-banner safe";
      alertBanner.innerHTML = `
        <div class="alert-icon">✅</div>
        <div class="alert-content">
          <h4>SAFE: CONCENTRATION WITHIN ACCEPTABLE LIMITS</h4>
          <p>Concentration at monitoring station is <strong>${formatSci(p.c_ppm, 4)} ppm</strong>, safely below the safety limit of <strong>${inp.c_target_ppm.toFixed(2)} ppm</strong> (${p.hazard_ratio_percent.toFixed(1)}% of threshold).</p>
        </div>
      `;
    } else {
      alertBanner.className = "alert-banner zero";
      alertBanner.innerHTML = `
        <div class="alert-icon">⏳</div>
        <div class="alert-content">
          <h4>ZERO CONCENTRATION (AHEAD OF ADVECTION FRONT)</h4>
          <p>The pollutant cloud front is at <strong>${p.advection_front_m.toFixed(1)} m</strong>. The receptor at <strong>${inp.x_query.toFixed(1)} m</strong> will be reached at approximately <strong>t = ${p.t_arrival_sec.toFixed(1)} s</strong>.</p>
        </div>
      `;
    }
  }

  // 2. Metrics Cards
  const valPpm = document.getElementById("resPpm");
  const valMg = document.getElementById("resMgM3");
  const valKg = document.getElementById("resKgM3");
  const valStatus = document.getElementById("resStatusBadge");
  const valFront = document.getElementById("resFrontDist");
  const valArrival = document.getElementById("resArrival");
  const valSigmas = document.getElementById("resSigmas");
  const valPeak = document.getElementById("resPeak");

  if (valPpm) valPpm.innerText = formatSci(p.c_ppm, 4) + " ppm";
  if (valMg) valMg.innerText = formatSci(p.c_mg_m3, 4) + " mg/m³";
  if (valKg) valKg.innerText = p.c_kg_m3.toExponential(6) + " kg/m³";

  if (valStatus) {
    if (p.exceeds) {
      valStatus.innerText = "EXCEEDED";
      valStatus.className = "badge danger";
    } else {
      valStatus.innerText = "SAFE";
      valStatus.className = "badge success";
    }
  }

  if (valFront) {
    const label = inp.model_type === "plume" ? "Plume Front:" : "Puff Center:";
    valFront.innerText = `${label} ${p.advection_front_m.toFixed(1)} m`;
  }

  if (valArrival) {
    const diff = p.advection_front_m - inp.x_query;
    const text =
      diff >= 0
        ? `Arrived (Passed by ${(diff / inp.u_wind).toFixed(1)}s)`
        : `Arrival in ${Math.abs(diff / inp.u_wind).toFixed(1)}s (at t=${p.t_arrival_sec.toFixed(1)}s)`;
    valArrival.innerText = text;
  }

  if (valSigmas) {
    if (inp.model_type === "plume") {
      valSigmas.innerText = `σy = ${p.sigmas.sy.toFixed(2)} m | σz = ${p.sigmas.sz.toFixed(2)} m`;
    } else {
      valSigmas.innerText = `σx,y = ${p.sigmas.sy.toFixed(2)} m | σz = ${p.sigmas.sz.toFixed(2)} m`;
    }
  }

  if (valPeak) {
    valPeak.innerText = `Peak: ${formatSci(p.peak_ppm, 3)} ppm @ ${p.t_peak_sec.toFixed(1)}s`;
  }

  // 3. Render Active Visualization Tab
  renderActiveVisualization();

  // 4. Update Summary Report Text
  updateReportText(res);
}

function renderActiveVisualization() {
  if (!AppState.lastResult) return;
  if (AppState.activeTab === "contour") {
    render2DContour(AppState.lastResult);
  } else if (AppState.activeTab === "breakthrough") {
    renderBreakthroughChart(AppState.lastResult);
  } else if (AppState.activeTab === "profiles") {
    renderProfilesChart(AppState.lastResult);
  }
}

// ----------------------------------------------------------------------------
// VISUALIZATION 1: 2D SPATIAL HEATMAP / CONTOUR CANVAS
// ----------------------------------------------------------------------------

function render2DContour(res) {
  const canvas = document.getElementById("contourCanvas");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  // Clear canvas
  ctx.fillStyle = "#0c1322";
  ctx.fillRect(0, 0, width, height);

  const grid = res.grid;
  const inp = res.inputs;
  const pt = res.point;

  const xVec = grid.x_vec;
  const yVec = grid.y_vec;
  const cMat = grid.c_matrix;
  const maxPpm = grid.max_ppm;

  const xMin = xVec[0];
  const xMax = xVec[xVec.length - 1];
  const yMin = yVec[0];
  const yMax = yVec[yVec.length - 1];

  // Margins for axes
  const padL = 60;
  const padR = 85; // space for colorbar
  const padT = 30;
  const padB = 45;
  const plotW = width - padL - padR;
  const plotH = height - padT - padB;

  function toScreenX(x) {
    return padL + ((x - xMin) / (xMax - xMin)) * plotW;
  }
  function toScreenY(y) {
    // Invert y so positive crosswind is up
    return padT + ((yMax - y) / (yMax - yMin)) * plotH;
  }
  function toDataX(px) {
    return xMin + ((px - padL) / plotW) * (xMax - xMin);
  }
  function toDataY(py) {
    return yMax - ((py - padT) / plotH) * (yMax - yMin);
  }

  // 1. Draw Grid Cells with Bilinear Color Interpolation
  const nx = xVec.length;
  const ny = yVec.length;
  const cellW = plotW / (nx - 1);
  const cellH = plotH / (ny - 1);

  // Color map function (0 to 1 -> color)
  function getColor(valPpm) {
    if (valPpm <= 1e-9) return "rgba(12, 19, 34, 0.0)";
    // Normalize relative to display peak
    const norm = Math.min(Math.max(valPpm / Math.max(maxPpm, 1e-6), 0.0), 1.0);

    // Color gradient stops: Navy -> Blue -> Teal -> Green -> Yellow -> Orange -> Crimson
    let r = 0, g = 0, b = 0, a = 0.88;
    if (norm < 0.05) {
      // near-zero faint cyan
      r = 30; g = 80; b = 150; a = norm / 0.05 * 0.45;
    } else if (norm < 0.25) {
      const f = (norm - 0.05) / 0.20;
      r = Math.round(20 + 20 * f);
      g = Math.round(90 + 120 * f);
      b = Math.round(180 + 40 * f);
      a = 0.5 + 0.35 * f;
    } else if (norm < 0.55) {
      const f = (norm - 0.25) / 0.30;
      r = Math.round(40 + 180 * f);
      g = Math.round(210 + 35 * f);
      b = Math.round(220 - 180 * f);
      a = 0.85;
    } else if (norm < 0.8) {
      const f = (norm - 0.55) / 0.25;
      r = Math.round(220 + 35 * f);
      g = Math.round(245 - 110 * f);
      b = Math.round(40 - 20 * f);
      a = 0.92;
    } else {
      const f = (norm - 0.8) / 0.20;
      r = 255;
      g = Math.round(135 - 100 * f);
      b = Math.round(20 + 50 * f);
      a = 0.98;
    }
    return `rgba(${r}, ${g}, ${b}, ${a})`;
  }

  // Draw heatmap pixels/rectangles
  for (let iy = 0; iy < ny - 1; iy++) {
    const y0 = toScreenY(yVec[iy]);
    const y1 = toScreenY(yVec[iy + 1]);
    const h = Math.abs(y0 - y1);
    const topY = Math.min(y0, y1);

    for (let ix = 0; ix < nx - 1; ix++) {
      const x0 = toScreenX(xVec[ix]);
      const x1 = toScreenX(xVec[ix + 1]);
      const w = Math.abs(x1 - x0);

      // Average cell value
      const val = (cMat[iy][ix] + cMat[iy][ix + 1] + cMat[iy + 1][ix] + cMat[iy + 1][ix + 1]) * 0.25;
      if (val > 1e-7) {
        ctx.fillStyle = getColor(val);
        ctx.fillRect(x0 - 0.5, topY - 0.5, w + 1.0, h + 1.0);
      }
    }
  }

  // 2. Draw Hazard Threshold Contour Line (Marching-squares approximation)
  const thresh = inp.c_target_ppm;
  if (maxPpm >= thresh) {
    ctx.strokeStyle = "#ff2a5f";
    ctx.lineWidth = 2.4;
    ctx.setLineDash([5, 4]);

    for (let iy = 0; iy < ny - 1; iy++) {
      for (let ix = 0; ix < nx - 1; ix++) {
        const v00 = cMat[iy][ix];
        const v10 = cMat[iy][ix + 1];
        const v11 = cMat[iy + 1][ix + 1];
        const v01 = cMat[iy + 1][ix];

        const minV = Math.min(v00, v10, v11, v01);
        const maxV = Math.max(v00, v10, v11, v01);

        if (thresh >= minV && thresh <= maxV) {
          // Cell crosses the threshold
          const x0 = toScreenX(xVec[ix]);
          const x1 = toScreenX(xVec[ix + 1]);
          const y0 = toScreenY(yVec[iy]);
          const y1 = toScreenY(yVec[iy + 1]);

          ctx.beginPath();
          ctx.rect(x0, Math.min(y0, y1), x1 - x0, Math.abs(y1 - y0));
          ctx.stroke();
        }
      }
    }
    ctx.setLineDash([]);
  }

  // 3. Gridlines & Axes
  ctx.strokeStyle = "rgba(255, 255, 255, 0.12)";
  ctx.lineWidth = 1;

  // Horizontal centerline (y = 0)
  const centerScreenY = toScreenY(0);
  ctx.beginPath();
  ctx.moveTo(padL, centerScreenY);
  ctx.lineTo(padL + plotW, centerScreenY);
  ctx.strokeStyle = "rgba(255, 255, 255, 0.25)";
  ctx.setLineDash([4, 4]);
  ctx.stroke();
  ctx.setLineDash([]);

  // Plot border
  ctx.strokeStyle = "rgba(255, 255, 255, 0.25)";
  ctx.strokeRect(padL, padT, plotW, plotH);

  // X-axis ticks & labels
  ctx.fillStyle = "#9aa5c4";
  ctx.font = "11px Consolas, monospace";
  ctx.textAlign = "center";
  const numXTicks = 6;
  for (let i = 0; i <= numXTicks; i++) {
    const xVal = xMin + (i / numXTicks) * (xMax - xMin);
    const sx = toScreenX(xVal);
    ctx.beginPath();
    ctx.moveTo(sx, padT + plotH);
    ctx.lineTo(sx, padT + plotH + 5);
    ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
    ctx.stroke();
    ctx.fillText(`${Math.round(xVal)}m`, sx, padT + plotH + 18);
  }
  ctx.fillText("Downwind Distance x (meters)", padL + plotW / 2, padT + plotH + 36);

  // Y-axis ticks & labels
  ctx.textAlign = "right";
  const numYTicks = 5;
  for (let i = 0; i <= numYTicks; i++) {
    const yVal = yMin + (i / numYTicks) * (yMax - yMin);
    const sy = toScreenY(yVal);
    ctx.beginPath();
    ctx.moveTo(padL - 5, sy);
    ctx.lineTo(padL, sy);
    ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
    ctx.stroke();
    ctx.fillText(`${yVal.toFixed(0)}m`, padL - 8, sy + 4);
  }

  // Y-axis title (rotated)
  ctx.save();
  ctx.translate(16, padT + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText("Crosswind y (meters)", 0, 0);
  ctx.restore();

  // 4. Mark Advection Wavefront / Center
  const frontX = pt.advection_front_m;
  if (frontX >= xMin && frontX <= xMax) {
    const sFrontX = toScreenX(frontX);
    ctx.strokeStyle = "#4cc9f0";
    ctx.lineWidth = 2;
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    ctx.moveTo(sFrontX, padT);
    ctx.lineTo(sFrontX, padT + plotH);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = "#4cc9f0";
    ctx.font = "10px Segoe UI, sans-serif";
    ctx.textAlign = "center";
    const label = inp.model_type === "plume" ? `Plume Front (${frontX.toFixed(0)}m)` : `Puff Center (${frontX.toFixed(0)}m)`;
    ctx.fillText(label, sFrontX, padT - 8);
  }

  // 5. Mark Release Source (0, 0)
  const sourceSX = toScreenX(0);
  const sourceSY = toScreenY(0);
  if (sourceSX >= padL && sourceSX <= padL + plotW) {
    ctx.fillStyle = "#f72585";
    ctx.beginPath();
    ctx.arc(sourceSX, sourceSY, 6, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = "#ffffff";
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 10px Segoe UI, sans-serif";
    ctx.textAlign = "left";
    ctx.fillText(" Source (0,0)", sourceSX + 7, sourceSY + 4);
  }

  // 6. Mark Monitoring Query Station (x_query, y_query)
  const querySX = toScreenX(inp.x_query);
  const querySY = toScreenY(inp.y_query);
  if (querySX >= padL && querySX <= padL + plotW && querySY >= padT && querySY <= padT + plotH) {
    // Crosshair circle
    ctx.fillStyle = pt.exceeds ? "#ff2a5f" : "#06d6a0";
    ctx.beginPath();
    ctx.arc(querySX, querySY, 7, 0, 2 * Math.PI);
    ctx.fill();
    ctx.strokeStyle = "#ffffff";
    ctx.lineWidth = 2;
    ctx.stroke();

    // Pulse effect ring
    ctx.strokeStyle = pt.exceeds ? "rgba(255, 42, 95, 0.5)" : "rgba(6, 214, 160, 0.5)";
    ctx.beginPath();
    ctx.arc(querySX, querySY, 13, 0, 2 * Math.PI);
    ctx.stroke();

    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 11px Segoe UI, sans-serif";
    ctx.textAlign = inp.x_query > (xMax + xMin) / 2 ? "right" : "left";
    const offset = inp.x_query > (xMax + xMin) / 2 ? -15 : 15;
    ctx.fillText(
      `Station (${inp.x_query.toFixed(0)}m, ${inp.y_query.toFixed(0)}m): ${formatSci(pt.c_ppm, 3)} ppm`,
      querySX + offset,
      querySY - 12
    );
  }

  // 7. Draw Colorbar Legend
  const cbL = width - padR + 15;
  const cbW = 16;
  const cbT = padT;
  const cbH = plotH;

  const cbSteps = 60;
  for (let i = 0; i <= cbSteps; i++) {
    const val = (1.0 - i / cbSteps) * maxPpm;
    const sy = cbT + (i / cbSteps) * cbH;
    ctx.fillStyle = getColor(val);
    ctx.fillRect(cbL, sy, cbW, cbH / cbSteps + 1);
  }
  ctx.strokeStyle = "rgba(255,255,255,0.3)";
  ctx.strokeRect(cbL, cbT, cbW, cbH);

  // Colorbar labels
  ctx.fillStyle = "#e8ecf8";
  ctx.font = "10px Consolas, monospace";
  ctx.textAlign = "left";
  ctx.fillText(`${formatSci(maxPpm, 2)}`, cbL + cbW + 6, cbT + 10);
  ctx.fillText(`${formatSci(maxPpm * 0.5, 2)}`, cbL + cbW + 6, cbT + cbH / 2 + 4);
  ctx.fillText(`0.00 ppm`, cbL + cbW + 6, cbT + cbH);

  ctx.save();
  ctx.translate(cbL + cbW + 42, cbT + cbH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.font = "11px Segoe UI, sans-serif";
  ctx.fillText("ppm", 0, 0);
  ctx.restore();

  // 8. Title Header
  ctx.fillStyle = "#e8ecf8";
  ctx.font = "bold 13px Segoe UI, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(
    `2D Concentration Contour Map at t = ${inp.t_query.toFixed(1)} s (z = ${inp.z_query.toFixed(1)} m) — Peak: ${formatSci(maxPpm, 2)} ppm`,
    padL,
    18
  );

  // Store coordinates mapping on canvas element for interactive hover tooltip
  canvas._toDataX = toDataX;
  canvas._toDataY = toDataY;
  canvas._padL = padL;
  canvas._padR = padR;
  canvas._padT = padT;
  canvas._padB = padB;
}

// ----------------------------------------------------------------------------
// VISUALIZATION 2: BREAKTHROUGH TIME-SERIES C(t)
// ----------------------------------------------------------------------------

function renderBreakthroughChart(res) {
  const canvas = document.getElementById("timeCanvas");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  ctx.fillStyle = "#0c1322";
  ctx.fillRect(0, 0, width, height);

  const tVec = res.temporal.t_vec;
  const cVec = res.temporal.c_vec;
  const inp = res.inputs;
  const pt = res.point;
  const thresh = inp.c_target_ppm;

  const tMin = tVec[0];
  const tMax = tVec[tVec.length - 1];
  const maxC = Math.max(...cVec, thresh * 1.25, 0.01);

  const padL = 65;
  const padR = 30;
  const padT = 36;
  const padB = 45;
  const plotW = width - padL - padR;
  const plotH = height - padT - padB;

  function toX(t) {
    return padL + ((t - tMin) / (tMax - tMin)) * plotW;
  }
  function toY(c) {
    return padT + ((maxC - c) / maxC) * plotH;
  }

  // Gridlines
  ctx.strokeStyle = "rgba(255, 255, 255, 0.08)";
  ctx.lineWidth = 1;
  const nTicks = 6;
  for (let i = 0; i <= nTicks; i++) {
    const tVal = tMin + (i / nTicks) * (tMax - tMin);
    const sx = toX(tVal);
    ctx.beginPath();
    ctx.moveTo(sx, padT);
    ctx.lineTo(sx, padT + plotH);
    ctx.stroke();

    ctx.fillStyle = "#9aa5c4";
    ctx.font = "11px Consolas, monospace";
    ctx.textAlign = "center";
    ctx.fillText(`${tVal.toFixed(0)}s`, sx, padT + plotH + 18);
  }

  for (let i = 0; i <= 5; i++) {
    const cVal = (i / 5) * maxC;
    const sy = toY(cVal);
    ctx.beginPath();
    ctx.moveTo(padL, sy);
    ctx.lineTo(padL + plotW, sy);
    ctx.stroke();

    ctx.fillStyle = "#9aa5c4";
    ctx.font = "11px Consolas, monospace";
    ctx.textAlign = "right";
    ctx.fillText(`${formatSci(cVal, 2)}`, padL - 8, sy + 4);
  }

  // Axis lines
  ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
  ctx.strokeRect(padL, padT, plotW, plotH);

  // Axis titles
  ctx.fillStyle = "#9aa5c4";
  ctx.font = "12px Segoe UI, sans-serif";
  ctx.textAlign = "center";
  ctx.fillText("Time Elapsed Since Release t (seconds)", padL + plotW / 2, padT + plotH + 36);

  ctx.save();
  ctx.translate(16, padT + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Concentration (ppm)", 0, 0);
  ctx.restore();

  // Hazard threshold line
  if (thresh <= maxC) {
    const threshY = toY(thresh);
    ctx.strokeStyle = "#ff2a5f";
    ctx.lineWidth = 2;
    ctx.setLineDash([6, 5]);
    ctx.beginPath();
    ctx.moveTo(padL, threshY);
    ctx.lineTo(padL + plotW, threshY);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = "#ff2a5f";
    ctx.font = "bold 11px Segoe UI, sans-serif";
    ctx.textAlign = "right";
    ctx.fillText(`Hazard Threshold (${thresh.toFixed(2)} ppm)`, padL + plotW - 10, threshY - 6);
  }

  // Evaluated time t_query vertical marker
  const evalX = toX(inp.t_query);
  ctx.strokeStyle = "#4cc9f0";
  ctx.lineWidth = 2;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(evalX, padT);
  ctx.lineTo(evalX, padT + plotH);
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.fillStyle = "#4cc9f0";
  ctx.font = "11px Segoe UI, sans-serif";
  ctx.textAlign = "center";
  ctx.fillText(`t = ${inp.t_query.toFixed(1)}s`, evalX, padT - 8);

  // Fill area under curve
  ctx.beginPath();
  ctx.moveTo(toX(tVec[0]), toY(0));
  for (let i = 0; i < tVec.length; i++) {
    ctx.lineTo(toX(tVec[i]), toY(cVec[i]));
  }
  ctx.lineTo(toX(tVec[tVec.length - 1]), toY(0));
  ctx.closePath();
  const fillGrad = ctx.createLinearGradient(0, padT, 0, padT + plotH);
  fillGrad.addColorStop(0, "rgba(76, 201, 240, 0.45)");
  fillGrad.addColorStop(1, "rgba(76, 201, 240, 0.02)");
  ctx.fillStyle = fillGrad;
  ctx.fill();

  // Draw C(t) Curve
  ctx.strokeStyle = "#4cc9f0";
  ctx.lineWidth = 2.8;
  ctx.beginPath();
  for (let i = 0; i < tVec.length; i++) {
    const sx = toX(tVec[i]);
    const sy = toY(cVec[i]);
    if (i === 0) ctx.moveTo(sx, sy);
    else ctx.lineTo(sx, sy);
  }
  ctx.stroke();

  // Mark current query point (t_query, C)
  const currentSY = toY(pt.c_ppm);
  ctx.fillStyle = pt.exceeds ? "#ff2a5f" : "#06d6a0";
  ctx.beginPath();
  ctx.arc(evalX, currentSY, 6, 0, 2 * Math.PI);
  ctx.fill();
  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 2;
  ctx.stroke();

  // Title header
  ctx.fillStyle = "#e8ecf8";
  ctx.font = "bold 13px Segoe UI, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(
    `Breakthrough Curve C(t) at Station (x = ${inp.x_query.toFixed(0)}m, y = ${inp.y_query.toFixed(0)}m, z = ${inp.z_query.toFixed(1)}m)`,
    padL,
    18
  );
}

// ----------------------------------------------------------------------------
// VISUALIZATION 3: CROSSWIND & DOWNWIND PROFILES
// ----------------------------------------------------------------------------

function renderProfilesChart(res) {
  const canvas = document.getElementById("profileCanvas");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  ctx.fillStyle = "#0c1322";
  ctx.fillRect(0, 0, width, height);

  const inp = res.inputs;
  const pt = res.point;
  const isCrosswind = AppState.profileSubTab === "crosswind";

  const xData = isCrosswind ? res.crosswind.y_vec : res.downwind.x_vec;
  const yData = isCrosswind ? res.crosswind.c_vec : res.downwind.c_vec;
  const xMin = xData[0];
  const xMax = xData[xData.length - 1];
  const maxC = Math.max(...yData, inp.c_target_ppm * 1.15, 0.01);

  const padL = 65;
  const padR = 30;
  const padT = 36;
  const padB = 45;
  const plotW = width - padL - padR;
  const plotH = height - padT - padB;

  function toX(x) {
    return padL + ((x - xMin) / (xMax - xMin)) * plotW;
  }
  function toY(c) {
    return padT + ((maxC - c) / maxC) * plotH;
  }

  // Grid
  ctx.strokeStyle = "rgba(255, 255, 255, 0.08)";
  ctx.lineWidth = 1;
  const nTicks = 6;
  for (let i = 0; i <= nTicks; i++) {
    const v = xMin + (i / nTicks) * (xMax - xMin);
    const sx = toX(v);
    ctx.beginPath();
    ctx.moveTo(sx, padT);
    ctx.lineTo(sx, padT + plotH);
    ctx.stroke();

    ctx.fillStyle = "#9aa5c4";
    ctx.font = "11px Consolas, monospace";
    ctx.textAlign = "center";
    ctx.fillText(`${v.toFixed(0)}m`, sx, padT + plotH + 18);
  }

  for (let i = 0; i <= 5; i++) {
    const cVal = (i / 5) * maxC;
    const sy = toY(cVal);
    ctx.beginPath();
    ctx.moveTo(padL, sy);
    ctx.lineTo(padL + plotW, sy);
    ctx.stroke();

    ctx.fillStyle = "#9aa5c4";
    ctx.font = "11px Consolas, monospace";
    ctx.textAlign = "right";
    ctx.fillText(`${formatSci(cVal, 2)}`, padL - 8, sy + 4);
  }

  // Border
  ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
  ctx.strokeRect(padL, padT, plotW, plotH);

  // Labels
  ctx.fillStyle = "#9aa5c4";
  ctx.font = "12px Segoe UI, sans-serif";
  ctx.textAlign = "center";
  const xLabel = isCrosswind
    ? "Crosswind Distance y (meters) at x = " + inp.x_query.toFixed(0) + "m"
    : "Centerline Downwind Distance x (meters) at y = 0";
  ctx.fillText(xLabel, padL + plotW / 2, padT + plotH + 36);

  ctx.save();
  ctx.translate(16, padT + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Concentration (ppm)", 0, 0);
  ctx.restore();

  // Hazard line
  if (inp.c_target_ppm <= maxC) {
    const threshY = toY(inp.c_target_ppm);
    ctx.strokeStyle = "#ff2a5f";
    ctx.lineWidth = 1.8;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.moveTo(padL, threshY);
    ctx.lineTo(padL + plotW, threshY);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Curve
  ctx.beginPath();
  ctx.moveTo(toX(xData[0]), toY(0));
  for (let i = 0; i < xData.length; i++) {
    ctx.lineTo(toX(xData[i]), toY(yData[i]));
  }
  ctx.lineTo(toX(xData[xData.length - 1]), toY(0));
  ctx.closePath();
  const fillGrad = ctx.createLinearGradient(0, padT, 0, padT + plotH);
  fillGrad.addColorStop(0, isCrosswind ? "rgba(255, 209, 102, 0.4)" : "rgba(6, 214, 160, 0.4)");
  fillGrad.addColorStop(1, "rgba(255, 209, 102, 0.02)");
  ctx.fillStyle = fillGrad;
  ctx.fill();

  ctx.strokeStyle = isCrosswind ? "#ffd166" : "#06d6a0";
  ctx.lineWidth = 2.6;
  ctx.beginPath();
  for (let i = 0; i < xData.length; i++) {
    const sx = toX(xData[i]);
    const sy = toY(yData[i]);
    if (i === 0) ctx.moveTo(sx, sy);
    else ctx.lineTo(sx, sy);
  }
  ctx.stroke();

  // Mark receptor
  const markVal = isCrosswind ? inp.y_query : inp.x_query;
  const markSX = toX(markVal);
  const markSY = toY(pt.c_ppm);
  ctx.fillStyle = pt.exceeds ? "#ff2a5f" : "#06d6a0";
  ctx.beginPath();
  ctx.arc(markSX, markSY, 6, 0, 2 * Math.PI);
  ctx.fill();
  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 2;
  ctx.stroke();

  // Title
  ctx.fillStyle = "#e8ecf8";
  ctx.font = "bold 13px Segoe UI, sans-serif";
  ctx.textAlign = "left";
  const title = isCrosswind
    ? `Crosswind Gaussian Profile C(y) at x = ${inp.x_query.toFixed(0)}m, t = ${inp.t_query.toFixed(1)}s`
    : `Downwind Centerline Profile C(x) at y = 0, t = ${inp.t_query.toFixed(1)}s`;
  ctx.fillText(title, padL, 18);
}

// ----------------------------------------------------------------------------
// 6. SUMMARY REPORT GENERATION
// ----------------------------------------------------------------------------

function updateReportText(res) {
  const box = document.getElementById("reportText");
  if (!box) return;
  const p = res.point;
  const inp = res.inputs;

  const text = `===========================================================================
   POINT CONCENTRATION REPORT AT t = ${inp.t_query.toFixed(1)} SECONDS
===========================================================================
Chemical & Model         : ${inp.chem_name} (${inp.model_type.toUpperCase()} Model)
Molecular Weight (MW)    : ${inp.mw_g_mol.toFixed(2)} g/mol
Release Strength         : ${inp.source_strength} ${inp.model_type === "plume" ? "kg/s" : "kg"}
Atmospheric Wind Speed   : ${inp.u_wind.toFixed(2)} m/s
Stability & Terrain      : Class ${inp.stability} | ${inp.terrain.toUpperCase()} terrain
Effective Release Height : Hr = ${inp.hr_m.toFixed(2)} m
---------------------------------------------------------------------------
Target Receptor Location : x = ${inp.x_query.toFixed(1)} m, y = ${inp.y_query.toFixed(1)} m, z = ${inp.z_query.toFixed(1)} m
Advection Front / Center : x_front = ${p.advection_front_m.toFixed(1)} m (at t = ${inp.t_query.toFixed(1)} s)
Expected Arrival Time    : t_arrival = ${p.t_arrival_sec.toFixed(1)} s
Dispersion Coefficients  : σy = ${p.sigmas.sy.toFixed(2)} m, σz = ${p.sigmas.sz.toFixed(2)} m
---------------------------------------------------------------------------
Concentration (kg/m³)    : ${p.c_kg_m3.toExponential(6)} kg/m³
Concentration (mg/m³)    : ${p.c_mg_m3.toFixed(4)} mg/m³
Concentration (ppm)      : ${p.c_ppm.toFixed(4)} ppm
Hazard Limit Threshold   : ${inp.c_target_ppm.toFixed(4)} ppm
Hazard Status            : ${p.exceeds ? "EXCEEDS THRESHOLD! (DANGER)" : "Below Threshold (SAFE)"}
Peak Concentration @ Rec : ${p.peak_ppm.toFixed(4)} ppm at t = ${p.t_peak_sec.toFixed(1)} s
===========================================================================`;

  box.innerText = text;
}

// ============================================================================
// 7. PRESET LOADER
// ============================================================================

function loadPreset(presetName) {
  const chemSelect = document.getElementById("chemSelect");
  const chemName = document.getElementById("chemName");
  const mwInput = document.getElementById("mwInput");
  const sourceStrength = document.getElementById("sourceStrength");
  const sourceUnit = document.getElementById("sourceUnit");
  const uWind = document.getElementById("uWind");
  const uWindSlider = document.getElementById("uWindSlider");
  const hrInput = document.getElementById("hrInput");
  const hrSlider = document.getElementById("hrSlider");
  const thresholdInput = document.getElementById("thresholdInput");
  const xQuery = document.getElementById("xQuery");
  const yQuery = document.getElementById("yQuery");
  const zQuery = document.getElementById("zQuery");
  const tQuery = document.getElementById("tQuery");
  const tQuerySlider = document.getElementById("tQuerySlider");

  if (presetName === "problem1_st1") {
    // Problem 1: SO2 continuous stack 80 g/s, u=6 m/s, Hr=60m, rural, class D, MW=64.07, threshold=0.5 ppm
    // Station 1: x=500m, y=0m, z=0m, t=100s
    chemSelect.value = "Sulfur Dioxide";
    chemName.value = "Sulfur Dioxide";
    mwInput.value = "64.07";
    document.querySelector('input[name="modelType"][value="plume"]').checked = true;
    updateModelTypeUI("plume");
    sourceStrength.value = "80";
    if (sourceUnit) sourceUnit.value = "g/s";
    uWind.value = "6.0";
    if (uWindSlider) uWindSlider.value = "6.0";
    document.querySelector('input[name="stability"][value="D"]').checked = true;
    document.querySelector('input[name="terrain"][value="rural"]').checked = true;
    hrInput.value = "60.0";
    if (hrSlider) hrSlider.value = "60";
    thresholdInput.value = "0.5";
    xQuery.value = "500.0";
    yQuery.value = "0.0";
    zQuery.value = "0.0";
    tQuery.value = "100.0";
    if (tQuerySlider) tQuerySlider.value = "100";
  } else if (presetName === "problem1_st2") {
    // Problem 1: Station 2: x=500m, y=50m, z=0m, t=100s
    loadPreset("problem1_st1");
    yQuery.value = "50.0";
  } else if (presetName === "problem2_st1") {
    // Problem 2: Chlorine instantaneous puff 100 kg, Hr=30m, u=5 m/s, rural, class D, MW=70.9, threshold=3.0 ppm
    // Station 1: x=500m, y=0m, z=0m, t=100s
    chemSelect.value = "Chlorine";
    chemName.value = "Chlorine";
    mwInput.value = "70.9";
    document.querySelector('input[name="modelType"][value="puff"]').checked = true;
    updateModelTypeUI("puff");
    sourceStrength.value = "100";
    if (sourceUnit) sourceUnit.value = "kg";
    uWind.value = "5.0";
    if (uWindSlider) uWindSlider.value = "5.0";
    document.querySelector('input[name="stability"][value="D"]').checked = true;
    document.querySelector('input[name="terrain"][value="rural"]').checked = true;
    hrInput.value = "30.0";
    if (hrSlider) hrSlider.value = "30";
    thresholdInput.value = "3.0";
    xQuery.value = "500.0";
    yQuery.value = "0.0";
    zQuery.value = "0.0";
    tQuery.value = "100.0";
    if (tQuerySlider) tQuerySlider.value = "100";
  } else if (presetName === "problem2_st2") {
    // Problem 2: Station 2: x=500m, y=50m, z=0m, t=100s
    loadPreset("problem2_st1");
    yQuery.value = "50.0";
  } else if (presetName === "default_cli") {
    // dispersion_model.py CLI default: Chlorine puff 50 kg, u=3 m/s, Hr=1.5 m, class D, rural, x=300m, y=0m, z=0m, t=100s
    chemSelect.value = "Chlorine";
    chemName.value = "Chlorine";
    mwInput.value = "70.9";
    document.querySelector('input[name="modelType"][value="puff"]').checked = true;
    updateModelTypeUI("puff");
    sourceStrength.value = "50";
    if (sourceUnit) sourceUnit.value = "kg";
    uWind.value = "3.0";
    if (uWindSlider) uWindSlider.value = "3.0";
    document.querySelector('input[name="stability"][value="D"]').checked = true;
    document.querySelector('input[name="terrain"][value="rural"]').checked = true;
    hrInput.value = "1.5";
    if (hrSlider) hrSlider.value = "1.5";
    thresholdInput.value = "3.0";
    xQuery.value = "300.0";
    yQuery.value = "0.0";
    zQuery.value = "0.0";
    tQuery.value = "100.0";
    if (tQuerySlider) tQuerySlider.value = "100";
  }

  // Highlight active preset button
  document.querySelectorAll(".preset-btn").forEach((b) => b.classList.remove("active"));
  const btn = document.querySelector(`[data-preset="${presetName}"]`);
  if (btn) btn.classList.add("active");

  setStabilityMode("direct");
  updateDirectStabilityHint("D");

  runSimulation();
}

function updateModelTypeUI(type) {
  const strengthLabel = document.getElementById("strengthLabel");
  const sourceUnit = document.getElementById("sourceUnit");
  if (type === "plume") {
    if (strengthLabel) strengthLabel.innerText = "Continuous Release Rate (Qm):";
    if (sourceUnit) {
      sourceUnit.innerHTML = `
        <option value="kg/s" selected>kg/s</option>
        <option value="g/s">g/s</option>
      `;
    }
  } else {
    if (strengthLabel) strengthLabel.innerText = "Total Released Mass (Qm*):";
    if (sourceUnit) {
      sourceUnit.innerHTML = `
        <option value="kg" selected>kg</option>
        <option value="g">g</option>
      `;
    }
  }
}

// ============================================================================
// 8. EVENT LISTENERS & INITIALIZATION
// ============================================================================

document.addEventListener("DOMContentLoaded", () => {
  // Check backend server availability
  // checkBackend();

  // Chemical dropdown change
  const chemSelect = document.getElementById("chemSelect");
  const chemName = document.getElementById("chemName");
  const mwInput = document.getElementById("mwInput");
  const thresholdInput = document.getElementById("thresholdInput");

  if (chemSelect) {
    chemSelect.addEventListener("change", () => {
      const val = chemSelect.value;
      if (val !== "Custom" && CHEMICAL_PRESETS[val]) {
        chemName.value = val;
        mwInput.value = CHEMICAL_PRESETS[val].mw;
        thresholdInput.value = CHEMICAL_PRESETS[val].defaultThreshold;
      }
      if (AppState.autoCalculate) runSimulation();
    });
  }

  // Model Type switch
  document.querySelectorAll('input[name="modelType"]').forEach((r) => {
    r.addEventListener("change", (e) => {
      updateModelTypeUI(e.target.value);
      if (AppState.autoCalculate) runSimulation();
    });
  });

  // Synchronize Sliders and Inputs
  function syncInputs(numInputId, sliderId) {
    const num = document.getElementById(numInputId);
    const slider = document.getElementById(sliderId);
    if (!num || !slider) return;

    slider.addEventListener("input", () => {
      num.value = slider.value;
      if (numInputId === "uWind" && AppState.stabilityMode === "auto") {
        updateEnvironmentalPrediction();
      }
      if (AppState.autoCalculate) runSimulation();
    });
    num.addEventListener("input", () => {
      slider.value = num.value;
      if (numInputId === "uWind" && AppState.stabilityMode === "auto") {
        updateEnvironmentalPrediction();
      }
      if (AppState.autoCalculate) runSimulation();
    });
  }

  syncInputs("uWind", "uWindSlider");
  syncInputs("hrInput", "hrSlider");
  syncInputs("tQuery", "tQuerySlider");

  // Form input changes for instant recalculation
  const allInputs = [
    "chemName",
    "mwInput",
    "sourceStrength",
    "sourceUnit",
    "uWind",
    "hrInput",
    "thresholdInput",
    "xQuery",
    "yQuery",
    "zQuery",
    "tQuery",
  ];
  allInputs.forEach((id) => {
    const el = document.getElementById(id);
    if (el) {
      el.addEventListener("input", () => {
        if (id === "uWind" && AppState.stabilityMode === "auto") {
          updateEnvironmentalPrediction();
        }
        if (AppState.autoCalculate) runSimulation();
      });
      el.addEventListener("change", () => {
        if (id === "uWind" && AppState.stabilityMode === "auto") {
          updateEnvironmentalPrediction();
        }
        if (AppState.autoCalculate) runSimulation();
      });
    }
  });

  // Radio button changes (stability and terrain)
  document.querySelectorAll('input[name="stability"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      if (AppState.stabilityMode === "auto" && !radio._autoSetting) {
        // User explicitly picked a stability class: switch to direct mode
        setStabilityMode("direct");
      }
      updateDirectStabilityHint(radio.value);
      if (AppState.autoCalculate) runSimulation();
    });
  });

  document.querySelectorAll('input[name="terrain"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      if (AppState.autoCalculate) runSimulation();
    });
  });

  // Stability Determination Mode Switcher (Slider / Segmented Buttons)
  const modeToggle = document.getElementById("stabilityModeToggle");
  if (modeToggle) {
    modeToggle.addEventListener("change", (e) => {
      setStabilityMode(e.target.checked ? "auto" : "direct");
      if (AppState.autoCalculate) runSimulation();
    });
  }

  const btnModeDirect = document.getElementById("btnModeDirect");
  if (btnModeDirect) {
    btnModeDirect.addEventListener("click", () => {
      setStabilityMode("direct");
      if (AppState.autoCalculate) runSimulation();
    });
  }

  const btnModeAuto = document.getElementById("btnModeAuto");
  if (btnModeAuto) {
    btnModeAuto.addEventListener("click", () => {
      setStabilityMode("auto");
      if (AppState.autoCalculate) runSimulation();
    });
  }

  // Environmental condition radio changes (Auto-Predict Mode)
  document.querySelectorAll('input[name="envCondition"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      AppState.envCondition = radio.value;
      updateEnvironmentalPrediction();
      if (AppState.autoCalculate) runSimulation();
    });
  });

  // Toggle Table 5-1 Matrix Accordion
  const btnToggleMatrix = document.getElementById("btnToggleMatrix");
  const matrixContent = document.getElementById("matrixTableContent");
  if (btnToggleMatrix && matrixContent) {
    btnToggleMatrix.addEventListener("click", () => {
      const isVisible = matrixContent.style.display !== "none";
      matrixContent.style.display = isVisible ? "none" : "block";
      btnToggleMatrix.classList.toggle("open", !isVisible);
    });
  }

  // Interactive Click on Matrix Cells to auto-select that row & column
  document.querySelectorAll("#pgMatrixTable .td-cell").forEach((cell) => {
    cell.addEventListener("click", () => {
      const colKey = cell.getAttribute("data-col");
      const tr = cell.closest("tr");
      const rowIndex = parseInt(tr?.getAttribute("data-row"), 10);
      if (colKey && !isNaN(rowIndex)) {
        // Set env condition
        const condRadio = document.querySelector(`input[name="envCondition"][value="${colKey}"]`);
        if (condRadio) condRadio.checked = true;
        AppState.envCondition = colKey;

        // Set representative wind speed for that row
        const rowSpeeds = [1.5, 2.5, 3.5, 5.0, 7.5];
        const newSpeed = rowSpeeds[rowIndex];
        const uWindInput = document.getElementById("uWind");
        const uWindSlider = document.getElementById("uWindSlider");
        if (uWindInput) uWindInput.value = newSpeed.toFixed(1);
        if (uWindSlider) uWindSlider.value = newSpeed.toFixed(1);

        // Switch to auto mode and update prediction
        setStabilityMode("auto");
        updateEnvironmentalPrediction();
        if (AppState.autoCalculate) runSimulation();
      }
    });
  });

  // Initial stability hint
  updateDirectStabilityHint("D");

  // Preset buttons
  document.querySelectorAll("[data-preset]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const preset = btn.getAttribute("data-preset");
      loadPreset(preset);
      // Smooth scroll to simulator
      document.getElementById("simulator")?.scrollIntoView({ behavior: "smooth" });
    });
  });

  // Run Simulation Button
  const runBtn = document.getElementById("runBtn");
  if (runBtn) {
    runBtn.addEventListener("click", () => {
      runSimulation();
    });
  }

  // Auto-calc toggle
  const autoToggle = document.getElementById("autoCalcToggle");
  if (autoToggle) {
    autoToggle.addEventListener("change", (e) => {
      AppState.autoCalculate = e.target.checked;
      if (AppState.autoCalculate) runSimulation();
    });
  }

  // Visualization Tabs
  document.querySelectorAll(".viz-tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".viz-tab-btn").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      AppState.activeTab = btn.getAttribute("data-tab");

      const profileSubtabs = document.getElementById("profileSubtabs");
      if (profileSubtabs) {
        profileSubtabs.style.display = AppState.activeTab === "profiles" ? "flex" : "none";
      }

      document.querySelectorAll(".viz-panel").forEach((p) => p.classList.remove("active"));
      const panel = document.getElementById(`panel_${AppState.activeTab}`);
      if (panel) panel.classList.add("active");

      renderActiveVisualization();
    });
  });

  // Profile sub-tabs (crosswind vs downwind)
  document.querySelectorAll(".profile-subtab").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".profile-subtab").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      AppState.profileSubTab = btn.getAttribute("data-subtab");
      renderProfilesChart(AppState.lastResult);
    });
  });

  // Copy Report to Clipboard
  const copyBtn = document.getElementById("copyReportBtn");
  if (copyBtn) {
    copyBtn.addEventListener("click", () => {
      const report = document.getElementById("reportText")?.innerText || "";
      navigator.clipboard.writeText(report).then(() => {
        const originalText = copyBtn.innerText;
        copyBtn.innerText = "✓ Copied to Clipboard!";
        setTimeout(() => (copyBtn.innerText = originalText), 2200);
      });
    });
  }

  // Canvas Mouse Hover Tooltip on 2D Contour
  const contourCanvas = document.getElementById("contourCanvas");
  const tooltip = document.getElementById("contourTooltip");
  if (contourCanvas && tooltip) {
    contourCanvas.addEventListener("mousemove", (e) => {
      if (!contourCanvas._toDataX) return;
      const rect = contourCanvas.getBoundingClientRect();
      const scaleX = contourCanvas.width / rect.width;
      const scaleY = contourCanvas.height / rect.height;
      const px = (e.clientX - rect.left) * scaleX;
      const py = (e.clientY - rect.top) * scaleY;

      if (
        px < contourCanvas._padL ||
        px > contourCanvas.width - contourCanvas._padR ||
        py < contourCanvas._padT ||
        py > contourCanvas.height - contourCanvas._padB
      ) {
        tooltip.style.display = "none";
        return;
      }

      const dataX = contourCanvas._toDataX(px);
      const dataY = contourCanvas._toDataY(py);
      const inp = AppState.lastResult?.inputs;
      if (!inp) return;

      const c_kg = computeConcentration(
        inp.model_type,
        inp.source_strength,
        inp.u_wind,
        inp.hr_m,
        dataX,
        dataY,
        inp.z_query,
        inp.t_query,
        inp.stability,
        inp.terrain
      );
      const ppm = (c_kg * 24.45 * 1e6) / inp.mw_g_mol;

      tooltip.style.display = "block";
      tooltip.style.left = `${e.clientX + 14}px`;
      tooltip.style.top = `${e.clientY - 28}px`;
      tooltip.innerHTML = `
        <strong>Downwind x:</strong> ${dataX.toFixed(1)} m<br/>
        <strong>Crosswind y:</strong> ${dataY.toFixed(1)} m<br/>
        <strong>Concentration:</strong> ${formatSci(ppm, 3)} ppm
      `;
    });

    contourCanvas.addEventListener("mouseleave", () => {
      tooltip.style.display = "none";
    });
  }

  // Initial Run with Default CLI Preset
  loadPreset("default_cli");
});
