import { recipeAdapter } from "./recipe_adapter.mjs";
import { recipeStartingValues, recipeStartingEvent, calibrationPatch } from "./instrument_calibrations.mjs";

const byId = id => document.getElementById(id);

export class RecipeController {
  constructor({
    engine, state, audition, getRoutingController,
    buildControls, buildPageValues, onChanged,
  }) {
    this.engine = engine;
    this.state = state;
    this.audition = audition;
    this.getRoutingController = getRoutingController;
    this.buildControls = buildControls;
    this.buildPageValues = buildPageValues;
    this.onChanged = onChanged;
    this.savedStates = new Map();
  }

  populate() {
    const select = byId("instrument-recipe");
    select.replaceChildren(...this.engine.recipes.map(recipe => {
      const option = document.createElement("option");
      option.value = recipe.index;
      option.textContent = recipe.name;
      return option;
    }));
    select.value = this.state.recipeIndex;
    select.onchange = () => {
      this.remember();
      const recipeIndex = Number(select.value);
      const key = this.engine.recipes.find(
        item => item.index === recipeIndex)?.key;
      this.activate(recipeIndex, this.savedStates.get(key));
      this.onChanged();
    };
  }

  remember() {
    if (!this.state.patch) return;
    this.savedStates.set(this.state.recipeKey, {
      macros: [...this.state.macros],
      patch: recipeAdapter(this.state.recipeKey).withValues(
        this.state.patch, this.engine.parameters, this.state.macros),
      event: { ...this.state.event },
    });
  }

  activate(recipeIndex, saved = null) {
    this.engine.setRecipe(recipeIndex);
    const recipe = this.engine.recipes.find(
      entry => entry.index === recipeIndex);
    if (!recipe) throw new Error(`unknown recipe index ${recipeIndex}`);
    this.state.recipeIndex = recipeIndex;
    this.state.recipeKey = recipe.key;
    byId("instrument-recipe").value = recipeIndex;
    const values = saved?.macros ??
      recipeStartingValues(recipe.key, this.engine.parameters);
    this.state.macros.splice(0, this.state.macros.length, ...values);
    const startingPatch = recipeAdapter(recipe.key).create(this.engine.parameters, values);
    this.state.patch = structuredClone(saved?.patch ?? calibrationPatch(
      { parameter_preset: recipe.key === "drum.kick.v1" ? "kick" : null },
      this.engine.parameters, values, startingPatch));
    Object.assign(this.state.event, saved?.event ?? recipeStartingEvent(recipe.key) ?? {});
    this.buildControls(false);
    const routing = this.getRoutingController();
    routing?.setPatch(this.state.patch);
    routing?.refreshPresentation();
    this.buildPageValues();
    this.audition.setRecipe(
      recipeIndex, values, recipeAdapter(recipe.key).routing(this.state.patch),
    );
  }
}
