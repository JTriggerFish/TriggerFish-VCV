import { ModuleRoles, ModuleTypes } from "./percussion_registry.mjs";
import { recipeAdapter } from "./recipe_adapter.mjs";
import { RoutingView } from "./routing_view.mjs";

const byId = id => document.getElementById(id);

export class RoutingController {
  constructor({ state, engine, audition, scheduleRender, setStatus }) {
    this.state = state;
    this.engine = engine;
    this.audition = audition;
    this.scheduleRender = scheduleRender;
    this.setStatus = setStatus;
  }

  bind() {
    this.view = new RoutingView({
      compact: byId("routing-compact"),
      expanded: byId("routing-expanded"),
      dialog: byId("routing-dialog"),
      onLayoutChange: patch => { this.state.patch = patch; },
      onRouteToggle: connectionId => this.#toggleRoute(connectionId),
    });
    this.view.setPatch(this.state.patch);
    byId("routing-close").onclick = () => byId("routing-dialog").close();
    byId("routing-compact").onkeydown = event => {
      if (event.key !== "Enter" && event.key !== " ") return;
      event.preventDefault();
      byId("routing-dialog").showModal();
      this.view.renderExpanded();
    };
    this.refreshPresentation();
  }

  setPatch(patch) { this.view.setPatch(patch); }

  refreshPresentation() {
    byId("routing-recipe-label").textContent = this.engine.recipes.find(
      recipe => recipe.index === this.state.recipeIndex,
    )?.name ?? "Compiled recipe";
    const roles = new Set(this.state.patch.nodes.map(node =>
      ModuleTypes.get(node.type).role));
    byId("routing-legend").replaceChildren(...[...roles].map(role => {
      const item = document.createElement("span");
      const swatch = document.createElement("i");
      swatch.style.setProperty("--module-colour", ModuleRoles[role].colour);
      item.append(swatch, ModuleRoles[role].name);
      return item;
    }));
  }

  #toggleRoute(connectionId) {
    const connection = this.state.patch.connections.find(
      item => item.id === connectionId);
    if (!connection) return;
    const previous = connection.enabled !== false;
    connection.enabled = !previous;
    try {
      recipeAdapter(this.state.recipeKey).validate(this.state.patch);
    } catch (error) {
      connection.enabled = previous;
      this.setStatus(String(error));
      this.view.setPatch(this.state.patch);
      return;
    }
    this.view.setPatch(this.state.patch);
    this.audition.setRouting(
      recipeAdapter(this.state.recipeKey).routing(this.state.patch));
    this.scheduleRender(false);
  }
}
