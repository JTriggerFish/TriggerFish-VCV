// Errors persist independently of transient Rendering / Ready status messages.
let installed = false;
let count = 0;
export function reportError(error, context = "") {
  if (typeof document === "undefined") return;
  const banner = document.getElementById("error-banner");
  if (!banner) return;
  const message = error instanceof Error ? error.message : String(error);
  document.getElementById("error-message").textContent = context ? `${context}: ${message}` : message;
  document.getElementById("error-count").textContent = ++count > 1 ? `${count} errors reported` : "Error";
  banner.hidden = false;
}

export function installErrorBanner() {
  if (installed) return;
  installed = true;
  document.getElementById("dismiss-error").onclick = () => {
    document.getElementById("error-banner").hidden = true; count = 0;
  };
  window.addEventListener("error", event => {
    reportError(event.error ?? event.message ?? "Resource failed to load", "Browser");
  });
  window.addEventListener("unhandledrejection", event => reportError(event.reason, "Async operation"));
}
