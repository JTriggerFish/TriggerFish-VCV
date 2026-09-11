import Fits from "./texture_trials.fit.json" with {type: "json"};

export function mountTextureTrials(select, restore, status) {
  select.replaceChildren(new Option("Choose a gong texture trial…", ""),
    ...Fits.map(fit => new Option(fit.name, fit.id)));
  select.disabled = Fits.length === 0;
  select.onchange = async () => {
    const fit = Fits.find(item => item.id === select.value);
    if (!fit) return;
    select.disabled = true;
    try {
      await restore({fit});
      document.getElementById("snapshot-name").value = fit.name;
    } catch (error) {
      select.value = "";
      status(error);
    } finally {
      select.disabled = false;
    }
  };
}
