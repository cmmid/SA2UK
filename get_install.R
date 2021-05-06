
.args <- if (interactive()) c(
  "rpack.txt", "covidm", ".install"
) else commandArgs(trailingOnly = TRUE)

pcks <- readLines(.args[1])
need <- !sapply(pcks, require, character.only = TRUE)
if (length(pcks[need])) install.packages(pcks[need])
if (dim(installed.packages()[pcks,])[1] == length(pcks)) {
  writeLines("COMPLETE\n", tail(.args, 1))
} else stop("Unable to install all packages.")

cm_path <- .args[2]

#' force build of covidm
#' all subsequent uses should have cm_force_shared = T
cm_force_rebuild = T;
cm_build_verbose = T;
cm_force_shared = F;
cm_version = 2;
source(file.path(cm_path, "R", "covidm.R"))

fIs_observer <- "
if (x.size() && t >= x[4]) { // if past the intervention start
  if (t < P.changes.ch.front().values.size()-2) { // if there is a tomorrow
    // the fI *reduction* model
    auto fIs_reduction = [&](double obscases, double case0, double k, double baseline) {
      return (1-baseline)/(1+exp(-k*(obscases-case0))) + baseline;
    };
    auto ocases = dyn(\"cases\", t, {}, {})*x[3]; // Ia (prevalence of active infections) * ascertainment
    auto newmul = 1.0 - fIs_reduction(ocases, x[0], x[1], x[2]);
    P.changes.ch.front().values[t+1].assign(16, newmul);

    // fIa, fIp
    // P.changes.ch[0].values[t+1].assign(16, newmul);
    // P.changes.ch[1].values[t+1].assign(16, newmul);
    // fIs
    // P.changes.ch[2].values[t+1].assign(16, newmul+(1-newmul)*newmul);
    // adjust the fIs multiplier for tomorrow based on todays observed prevalence
  };
};
"

cm_source_backend(
  user_defined = list(
    model_v2 = list(
      cpp_changes = "",
      cpp_loglikelihood = "",
      cpp_observer = fIs_observer
    )
  )
)