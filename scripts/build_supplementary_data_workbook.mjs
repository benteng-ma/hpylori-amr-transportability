import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = process.env.HPYLORI_PROJECT_ROOT
  ? path.resolve(process.env.HPYLORI_PROJECT_ROOT)
  : path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const threadId = "01a05070-a483-7892-bc5a-ef5674322155";
const outputDir = process.env.HPYLORI_WORKBOOK_OUTPUT_DIR
  ? path.resolve(process.env.HPYLORI_WORKBOOK_OUTPUT_DIR)
  : path.join(root, "outputs", threadId);
const previewDir = path.join(outputDir, "workbook_previews");
const submissionPath = path.join(root, "submission", "Supplementary_Data_S1-S20.xlsx");
const tableDir = path.join(root, "results", "supplementary_tables");

const tables = [
  ["S01_Crosswalk", "S01_Crosswalk.csv"],
  ["S02_QC_Callability", "S02_QC_Callability.csv"],
  ["S03_Lineage_Relatedness", "S03_Lineage_Relatedness.csv"],
  ["S04_Predictions", "S04_Predictions.csv"],
  ["S05_Performance", "S05_Performance.csv"],
  ["S06_Validation_Design", "S06_Validation_Design.csv"],
  ["S07_Phenotype_Sensitivity", "S07_Phenotype_Sensitivity.csv"],
  ["S08_Negative_Controls", "S08_Negative_Controls.csv"],
  ["S09_MIC_Summary", "S09_MIC_Summary.csv"],
  ["S10_Genotype_Audit", "S10_Genotype_Audit.csv"],
  ["S11_Classification", "S11_Classification.csv"],
  ["S12_Figure_Manifest", "S12_Figure_Manifest.csv"],
  ["S13_Marker_Architecture", "S13_Marker_Architecture.csv"],
  ["S14_QC_Associations", "S14_QC_Associations.csv"],
  ["S15_Population_Shift", "S15_Population_Shift.csv"],
  ["S16_Residual_Mechanisms", "S16_Residual_Mechanisms.csv"],
  ["S17_Stability_Scenarios", "S17_Stability_Scenarios.csv"],
  ["S18_Transport_Shift", "S18_Transport_Shift.csv"],
  ["S19_23S_Sensitivity", "S19_23S_Sensitivity.csv"],
  ["S20_Coverage_Selective", "S20_Coverage_Selective.csv"],
];

function parseHeader(csvText) {
  const line = csvText.split(/\r?\n/, 1)[0];
  const fields = [];
  let current = "";
  let quoted = false;
  for (let index = 0; index < line.length; index += 1) {
    const character = line[index];
    if (character === '"') quoted = !quoted;
    else if (character === "," && !quoted) {
      fields.push(current.replace(/^"|"$/g, ""));
      current = "";
    } else current += character;
  }
  fields.push(current.replace(/^"|"$/g, ""));
  return fields;
}

function preferredWidth(header) {
  const lower = header.toLowerCase();
  if (lower.includes("source_url")) return 42;
  if (lower.includes("path") || lower.includes("source_file")) return 38;
  if (lower.includes("notes") || lower.includes("reason") || lower.includes("method")) return 34;
  if (lower.includes("checksum") || lower.includes("sha256")) return 24;
  if (lower.includes("status") || lower.includes("classification")) return 22;
  if (lower.includes("isolate") || lower.includes("dataset") || lower.includes("accession")) return 19;
  return 14;
}

await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.mergeCells("A1:H1");
readme.getRange("A1:H1").values = [["Supplementary Data S1-S20 — H. pylori AMR transportability benchmark"]];
readme.getRange("A1:H1").format = {
  fill: "#17365D",
  font: { bold: true, color: "#FFFFFF", size: 16 },
  verticalAlignment: "center",
};
readme.getRange("A1:H1").format.rowHeight = 34;
readme.getRange("A3:B13").values = [
  ["Purpose", "Consolidated, reviewer-ready machine-readable supplement for the audit-first transportability benchmark."],
  ["Scope", "526 materialized genomes; 414 HpGP assemblies, 60 Ningxia assemblies, and 52 Zenodo paired-read isolates."],
  ["Primary unit", "One row per isolate in S01; other sheets retain a source_file column when related outputs were stacked."],
  ["Frozen analyses", "S04-S11 contain the prespecified predictions, performance estimates, controls, and transportability labels."],
  ["Post-freeze analyses", "S13-S20 are explanatory or robustness analyses and do not alter frozen catalogues or classifications."],
  ["23S sensitivity", "S19 varies BLAST task, identity, and coverage while retaining the requirement that both resistance-marker bases are spanned."],
  ["Coverage/selective stress", "S20 reports end-to-end correct-result yield, logical missing-result bounds, phenotype-blind abstention, and external-cohort differences."],
  ["Missing values", "Blank cells mean unavailable, inapplicable, or not estimable; they are not recoded as susceptible or zero."],
  ["Provenance", "Each stacked record retains source_file. Stable raw-accession URLs and checksums are provided in S01."],
  ["Generated", "2026-09-01 (Asia/Shanghai project date)."],
  ["Authorship", "Five-author lock: Benteng Ma and Bing Chen contributed equally, with Benteng Ma listed first; Ting Cai is the only other author; Xiao-ming Liu and Fen Wang are co-corresponding authors."],
];
readme.getRange("A3:A13").format = { fill: "#D9EAF7", font: { bold: true, color: "#17365D" } };
readme.getRange("A3:B13").format.borders = { preset: "all", style: "thin", color: "#B7C9D6" };
readme.getRange("A3:B13").format.wrapText = true;
readme.getRange("A3:A13").format.columnWidth = 23;
readme.getRange("B3:B13").format.columnWidth = 95;
readme.getRange("A15:H15").merge();
readme.getRange("A15:H15").values = [["Worksheet index"]];
readme.getRange("A15:H15").format = { fill: "#2F75B5", font: { bold: true, color: "#FFFFFF" } };
readme.getRange("A16:C36").values = [
  ["Sheet", "Contents", "Interpretive status"],
  ["S01_Crosswalk", "Isolate–sequence–phenotype–prediction crosswalk", "Audit/provenance"],
  ["S02_QC_Callability", "Genome quality and target callability", "Frozen"],
  ["S03_Lineage_Relatedness", "Near clones, pairwise checks, fixed SNP clusters", "Frozen"],
  ["S04_Predictions", "Per-isolate frozen-panel predictions", "Frozen"],
  ["S05_Performance", "Cohort/country/lineage performance", "Frozen"],
  ["S06_Validation_Design", "Leakage-aware validation benchmark", "Frozen"],
  ["S07_Phenotype_Sensitivity", "Phenotype and breakpoint scenarios", "Sensitivity"],
  ["S08_Negative_Controls", "Clone thinning, concordance, random panels, permutations", "Negative controls"],
  ["S09_MIC_Summary", "Ningxia MIC lower-bound summaries", "Descriptive"],
  ["S10_Genotype_Audit", "Source-genotype and independent-caller audit", "Audit"],
  ["S11_Classification", "Prespecified transportability classification", "Frozen"],
  ["S12_Figure_Manifest", "Main-figure panel source-data map", "Provenance"],
  ["S13_Marker_Architecture", "Mutation spectra and prevalence shifts", "Post-freeze"],
  ["S14_QC_Associations", "QC–callability associations", "Post-freeze"],
  ["S15_Population_Shift", "Development-manifold and error-distance analyses", "Post-freeze"],
  ["S16_Residual_Mechanisms", "Ningxia residual-mechanism audit", "Hypothesis-generating"],
  ["S17_Stability_Scenarios", "Bootstrap, influence, prevalence and gate grids", "Robustness"],
  ["S18_Transport_Shift", "Four-domain transport-shift decomposition", "Post-freeze"],
  ["S19_23S_Sensitivity", "Target-recovery threshold/task audit", "Post-freeze robustness"],
  ["S20_Coverage_Selective", "End-to-end yield, logical bounds, abstention and external differences", "Post-freeze robustness"],
];
readme.getRange("A16:C16").format = { fill: "#4472C4", font: { bold: true, color: "#FFFFFF" } };
readme.getRange("A16:C36").format.borders = { preset: "all", style: "thin", color: "#D9E2F3" };
readme.getRange("A16:C36").format.wrapText = true;
readme.getRange("A16:A36").format.columnWidth = 27;
readme.getRange("B16:B36").format.columnWidth = 62;
readme.getRange("C16:C36").format.columnWidth = 25;
readme.freezePanes.freezeRows(2);

for (let index = 0; index < tables.length; index += 1) {
  const [sheetName, fileName] = tables[index];
  const csvText = await fs.readFile(path.join(tableDir, fileName), "utf8");
  const importSafeCsvText = csvText.replace(/(^|,)=(?=,|\r?$)/gm, "$1'=");
  const importedWorkbook = await Workbook.fromCSV(importSafeCsvText, { sheetName });
  const importedSheet = importedWorkbook.worksheets.getItem(sheetName);
  const importedValues = importedSheet.getUsedRange().values;
  const safeValues = importedValues.map((row) =>
    row.map((value) => (typeof value === "string" && value.startsWith("=") ? `'${value}` : value)),
  );
  const sheet = workbook.worksheets.add(sheetName);
  sheet.getRangeByIndexes(0, 0, safeValues.length, safeValues[0].length).values = safeValues;
  const used = sheet.getUsedRange();
  const headers = parseHeader(csvText);
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  used.format.font = { name: "Aptos", size: 9, color: "#1F1F1F" };
  used.format.verticalAlignment = "top";
  sheet.getRangeByIndexes(0, 0, 1, headers.length).format = {
    fill: index % 2 === 0 ? "#1F4E78" : "#2F5597",
    font: { bold: true, color: "#FFFFFF", size: 9 },
    wrapText: true,
    verticalAlignment: "center",
  };
  sheet.getRangeByIndexes(0, 0, 1, headers.length).format.rowHeight = 32;
  headers.forEach((header, column) => {
    used.getColumn(column).format.columnWidth = preferredWidth(header);
  });
  used.format.borders = { preset: "inside", style: "thin", color: "#E7E6E6" };
  sheet.tables.add(used, true, `T${String(index + 1).padStart(2, "0")}`);
}

const inspect = await workbook.inspect({
  kind: "workbook,sheet,table",
  maxChars: 12000,
  tableMaxRows: 3,
  tableMaxCols: 8,
  tableMaxCellChars: 60,
});
await fs.writeFile(path.join(outputDir, "Supplementary_Data_S1-S20.inspect.txt"), inspect.ndjson ?? String(inspect), "utf8");

const s20Inspect = await workbook.inspect({
  kind: "table",
  sheetId: "S20_Coverage_Selective",
  range: "A1:L30",
  include: "values,formulas",
  tableMaxRows: 30,
  tableMaxCols: 12,
  tableMaxCellChars: 120,
});
await fs.writeFile(
  path.join(outputDir, "Supplementary_Data_S1-S20_S20.inspect.txt"),
  s20Inspect.ndjson ?? String(s20Inspect),
  "utf8",
);

const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "Final formula-error scan",
});
await fs.writeFile(
  path.join(outputDir, "Supplementary_Data_S1-S20_formula_errors.inspect.txt"),
  formulaErrors.ndjson ?? String(formulaErrors),
  "utf8",
);

const exported = await SpreadsheetFile.exportXlsx(workbook);
const canonicalPath = path.join(outputDir, "Supplementary_Data_S1-S20.xlsx");
await exported.save(canonicalPath);
await fs.copyFile(canonicalPath, submissionPath);

for (const sheetName of ["README", ...tables.map(([name]) => name)]) {
  const preview = await workbook.render({ sheetName, range: "A1:L30", autoCrop: "all", scale: 0.8, format: "png" });
  await fs.writeFile(path.join(previewDir, `${sheetName}.png`), new Uint8Array(await preview.arrayBuffer()));
}
