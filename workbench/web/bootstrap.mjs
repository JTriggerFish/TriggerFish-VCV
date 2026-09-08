import {installErrorBanner, reportError} from "./error_banner.mjs";
installErrorBanner();
// Dynamic import lets the banner also report application module-load failures.
import("./app.mjs").catch(error => reportError(error, "Workbench startup"));
