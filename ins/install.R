
.args <- if (interactive()) c(
  "ins/rpack.csv", "covidm", ".install"
) else commandArgs(trailingOnly = TRUE)

pcks <- read.csv(.args[1], header = FALSE, strip.white = TRUE, col.names = c("source", "name"))
cranset <- subset(pcks, source == "cran")
install.packages(cranset$name, dependencies = TRUE, repos = "https://cloud.r-project.org/")

githubset <- subset(pcks, source == "github")
for (i in 1:nrow(githubset)) remotes::install_github(githubset$name[i])

cm_path <- .args[2]

#' force build of covidm
#' all subsequent uses should have cm_force_shared = T
cm_force_rebuild = T;
cm_build_verbose = T;
cm_force_shared = F;
cm_version = 2;
source(file.path(cm_path, "R", "covidm.R"))

fIs_observer <- "
double case0 = x[0];
double k = x[1];
double asc = x[2];
double modstartt = x[3];
double refsympt = x[4];
if (x.size() && t >= modstartt) { // if past the intervention start
  if (t < P.changes.ch.front().values.size()-2) { // if there is a tomorrow
    auto fIs_reduction = [&](double cases) {
      return 1/(1+exp(-k*(cases-case0)));
    };
    auto refred = fIs_reduction(dyn(\"cases\", modstartt, {}, {}));

    // the fI *reduction* model
    auto nowred = fIs_reduction(dyn(\"cases\", t, {}, {}));
    double newmul = 1.0;
    if (refsympt <= refred) {
      newmul -= refsympt/refred * nowred;
    } else {
      auto flr = (refsympt-refred)/(1-refred);
      newmul -= flr + (1-flr)*nowred;
    }
    
    P.changes.ch.front().values[t+1].assign(16, newmul);
//    std::cout << t << ' ' << ocases << ' ' << newmul << std::endl;
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
