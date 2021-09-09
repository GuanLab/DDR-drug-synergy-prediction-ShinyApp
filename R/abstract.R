#' @title 
#' Abstract
#' 
#' @description 
#' Abstract text for the landing page
#' 
abstract = function() {
  shiny::HTML("One principle underlying many types of cancer treatments,",
  "including traditional chemotherapy and radiotherapy, is inducing",
  "catastrophic DNA damage that results in apoptosis of cancer cells.",
  "Chemotherapy agents, as well as inhibitors that directly target DNA damage",
  "repair (DDR) pathways, show particularly strong efficacy in tumors that are",
  "vulnerable due to absence of certain DDR pathways, a concept known as",
  "synthetic lethality. In this study, we produced a comprehensive",
  "dose-response combination screen focusing on the DDR kinases",
  "Ataxia-telangiectasia mutated (ATM), Ataxia-telangiectasia and Rad3-related",
  "(ATR), and DNA-dependent protein kinase (DNA-PK), that are currently",
  "targeted in clinical development. The screen was conducted with 90 DDR and",
  "anti-cancer agents covering 62 cell lines across 12 cancer indications, thus",
  "forming a total of 17,912 combination treatment experiments. Analysis of",
  "this screen identified inhibitors of five DDR-related pathways (the DNA",
  "topoisomerase pathway, the serine/threonine-protein kinase PLK1 pathway,",
  "the p53-inducible ribonucleotide reductase pathway, the PARP pathway, and",
  "the cell cycle checkpoint proteins) that displayed particularly high",
  "combinatorial efficacy with ATM/ATR/DNA-PK inhibitors. Secondly, we found",
  "that the screened ATM, ATR, and DNA-PK inhibitors achieve strong synergistic",
  "effects with a PARP inhibitor that is used in first-line clinical cancer",
  "therapy today. Last, correlating dose responses with molecular readouts of",
  "these cell lines allowed us to develop predictive machine learning models",
  "integrating multi-scale, heterogeneous genomic and chemical information that",
  "improve treatment efficacy in over 95% of cases compared to baseline. Using",
  "game-theory-based feature analysis, we identified novel candidate molecular",
  "biomarkers that are predictive for cancer cell line models' sensitivity to",
  "DDR combination therapies.")
}
