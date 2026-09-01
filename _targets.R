library(targets)
library(jsonlite)

decision_path <- "reports/phase0_decision.json"
phase0_authorized <- function(path) {
  if (!file.exists(path)) return(FALSE)
  decision <- fromJSON(path)$decision
  grepl("^(GO_|CONDITIONAL_GO_)", decision)
}

tar_option_set(packages = c("data.table", "jsonlite", "yaml"))

list(
  tar_target(phase0_decision_path, decision_path, format = "file"),
  tar_target(phase0_gate, {
    if (!phase0_authorized(phase0_decision_path)) {
      stop("Full analysis is locked until Phase 0 records GO or CONDITIONAL_GO.")
    }
    TRUE
  })
)

