// TODO: pictures of galaxies as backgrounds?
// TODO: fix slightly odd behavior at top end of BH slider


const G    = 6.6743e-11; // N m^2/kg^2
const ly   = 9.4607e15; // m
const kpc  = 3.086e+19; // m
const MSun = 1.9891e30; // kg

const rUnit = 1000 * ly;
const mUnit = 1e9 * MSun;
const vUnit = 1000; // m/s

const alpha = 2.3;

const nPoints = 100;  // Number of points for plotting predicted velocities.


// List of objects representing galaxies, with one entry per galaxy in the
// selection menu.  Each entry is a simple object with the following keys:
//   name: the name of the galaxy
//   galaxy: its Galaxy object
//   model: its Model object
//
// This list is initialized by calling setupGalaxies().
var galaxyList;

// The Galaxy and Model instances for the currently selected galaxy.  These
// variables are references into galaxyList, so their changes persist if the
// user changes back and forth between selected galaxies.
//
// These are initialized by calling setupGalaxies().
var galaxy, model;

// Variables containing page elements.  These are initialized by onload().
var galaxySelector, velocityPlot, densityPlot;
var summaryText, bhSlider, bhText, showLMCheckbox, showLMLabel;


class Galaxy {
  // Input: g -- galaxy data object
  constructor(g) {
    this.name = g.name;
    this.rMin = g.rMin / rUnit;
    this.rMax = g.rMax / rUnit;
    this.rStep = g.rStep / rUnit;
    this.vMax = g.vMax / vUnit;
    this.mMax = g.mMax / mUnit;
    this.specialPoints =
        (g.specialPoints !== undefined) ? g.specialPoints : [];

    this.luminousMassModel = g.luminousMassModel;

    // Read in the observed velocities.
    this.r = [];
    this.v = [];
    this.vErr = [];
    this.minDataSpacing = 9e99;
    var rPrev = -9e99;
    for (var i = 0; i < g.data.length; i++) {
      var r = g.data[i][0] * g.rUnit / rUnit;
      var v = g.data[i][1] * g.vUnit / vUnit;
      var vErr = g.data[i][2] * g.vUnit / vUnit;
      this.r.push(r);
      this.v.push(v);
      this.vErr.push(vErr);
      if (this.minDataSpacing > r - rPrev) {
        this.minDataSpacing = r - rPrev;
      }
      rPrev = r;
    }
  }

  // Return an array for plotting this model on an RGraph scatterplot.
  makePlotData() {
    var data = [];
    for (var i = 0; i < this.r.length; i++) {
      var v = this.v[i], vErr = this.vErr[i];
      data.push([this.r[i], [v - vErr, v, v, v, v + vErr]]);
    }
    return data;
  }

  getLuminousMassCurve(radii) {
    if (this.luminousMassModel === undefined) {
      return [];
    }

    const m0 = this.luminousMassModel.m0;
    const r0 = this.luminousMassModel.r0;
    const mu = this.luminousMassModel.mu;
    const nu = this.luminousMassModel.nu;
    var shellMasses = [];
    for (var r of radii) {
      shellMasses.push(
          m0 * Math.pow(r / r0, mu) * Math.exp(-Math.pow(r / r0, nu)));
    }
    return shellMasses;
  }
}


// Our model:
// * Assume mass is in shells from r_{i-1} to r_i.
// * Shell i = 0 is a central point mass.
// * All other shells are solid with mass accummulating as r^alpha.
// * For constant density spherical shells, r = 3; for disk, r = 2.
// * We choose to use r = 2.3 because it produces nice-looking plots.
// * r_i goes from 0 to rMax in steps of rStep.

class Model {
  // Given a Galaxy object, initialize:
  // * a mass model, set to zero.
  // * an array of radii for plotting a predicted veloctiy curve.
  //
  // Input: g -- Galaxy object to be modeled.
  constructor(g) {
    this.r = [];        // Radii at which our model is defined, in rUnit.
    this.rLabels = [];  // Text labels for the above radii.
    this.m = [];        // Masses for shells ending at these radii, in mUnit.
    for (var r = 0; r < g.rMax + 0.1 * g.rStep; r += g.rStep) {
      this.r.push(r);
      this.rLabels.push(r.toString())
      this.m.push(0);
    }
    this.rPlot = [];
    for (var i = 0; i < nPoints; i++) {
      this.rPlot.push((g.rMax - g.rMin) * i / (nPoints - 1) + g.rMin);
    }
  }

  // Update the model based on the adjustable plot data.
  //
  // Input: data -- array of masses in units of mUnit.
  update(data) {
    var mTotal = 0;
    for (var i = 0; i < this.m.length; i++) {
      this.m[i] = data[i];
      mTotal += data[i];
    }
    return mTotal;
  }

  // Use this model to predict velocities at the given radii.
  //
  // Input: radii -- an array of radii, in units of rUnit.
  //
  // Returns: an array of velocities, in units of vUnit.
  calcVelocities(radii) {
    var vels = [], i = 0, rPrev = 0, mPrev = 0;
    for (var r of radii) {
      while (this.r[i] <= r) {
        mPrev += this.m[i];
        rPrev = this.r[i];
        i++;
      }
      var m = mPrev;
      if (i < this.m.length) {
        m += this.m[i] * (
            (Math.pow(r, alpha) - Math.pow(rPrev, alpha)) /
                (Math.pow(this.r[i], alpha) - Math.pow(rPrev, alpha)));
      }
      vels.push(Math.sqrt(G * (m * mUnit) / (r * rUnit)) / vUnit);
    }
    return vels;
  }

  // Given an observed rotation curve, calculate the reduced chi squared for
  // this model.
  //
  // Input: galaxy -- Galaxy object that contains the observed rotation curve.
  //
  // Returns: reduced chi squared statistic.
  calcReducedChiSqrd(galaxy) {
    var nDoF = galaxy.r.length - this.r.length;
    if (nDoF <= 0) {
      return -1;
    }
    var vPred = this.calcVelocities(galaxy.r);
    var chiSqrd = 0;
    for (var i = 0; i < vPred.length; i++) {
      chiSqrd += Math.pow((galaxy.v[i] - vPred[i]) / galaxy.vErr[i], 2);
    }
    return chiSqrd / nDoF;
  }

  // Return an array for plotting this model on an RGraph scatterplot.
  makePlotData() {
    var vels = this.calcVelocities(this.rPlot);
    var data = [];
    for (var i = 0; i < this.rPlot.length; i++) {
      data.push([this.rPlot[i], vels[i]]);
    }
    return data;
  }
}


// Populate galaxyList by creating Galaxy and Model instances for the given
// galaxy data.  Set the initial galaxy to the first one on the list.
function setupGalaxies(galaxyData) {
  galaxyList = [];
  for (gd of galaxyData) {
    var g = new Galaxy(gd);
     galaxyList.push({
      name: gd.name,
      galaxy: g,
      model: new Model(g)
    });
  }
  galaxy = galaxyList[0].galaxy;
  model = galaxyList[0].model;
}


// Initialize DOM variables and setup the page.
window.onload = function() {
  galaxySelector = document.getElementById("galaxy_selector");
  summaryText = document.getElementById("summary");
  bhSlider = document.getElementById("bh_slider");
  bhText = document.getElementById("bh_mass");
  showLMCheckbox = document.getElementById("show_lum_mass")
  showLMLabel = document.getElementById("show_lum_mass_label")

  for (var i = 0; i < galaxyList.length; i++) {
    var opt = document.createElement("option");
    opt.value = i;
    opt.innerHTML = galaxyList[i].name;
    galaxySelector.appendChild(opt);
  }
  galaxySelector.value = 0;
  galaxySelector.onchange = onUpdateGalaxySelection;

  bhSlider.oninput = onUpdateBhSlider;
  showLMCheckbox.onclick = onShowLumMassClicked;

  setupPlots();
}


function onUpdateGalaxySelection() {
  var i = galaxySelector.value;
  galaxy = galaxyList[i].galaxy;
  model = galaxyList[i].model;
  setupPlots();
}


function setupPlots() {
  if (velocityPlot) {
    RGraph.Clear(velocityPlot.canvas);
    RGraph.Clear(densityPlot.canvas);
    RGraph.ObjectRegistry.Clear();
  }

  velocityPlot = new RGraph.Scatter({
    id: "velocity",
    data: [galaxy.makePlotData(), model.makePlotData()],
    options: {
      xaxisScaleMin: 0,
      xaxisScaleMax: galaxy.rMax,
      yaxisScaleMax: galaxy.vMax,
      yaxisTitle: "Orbital velocity (km/s)",
      yaxisTitlePos: 0.3,
      xaxisTickmarksCount: model.r.length - 1,
      backgroundGridVlinesCount: model.r.length - 1,
      marginLeft: 50,
      marginRight: 10,
      marginBottom: 5,
      line: true,
      boxplotWidth: 0.75 * galaxy.minDataSpacing,
      boxplotCapped: true,
      tickmarksStyle: [null, null, "circle"],
      tickmarksSize: 5,
      textFont: "Arial",
    }
  }).draw();

  for (var pt of galaxy.specialPoints) {
    new RGraph.Drawing.Circle({
      id: "velocity",
      x: velocityPlot.getXCoord(pt.r),
      y: velocityPlot.getYCoord(pt.v),
      radius: pt.ptRad,
      options: {
        colorsStroke: "black",
        colorsFill: pt.color,
      }
    }).draw();
  }

  densityPlot = new RGraph.Line({
    id: "density",
    data: [model.m, galaxy.getLuminousMassCurve(model.r)],
    options: {
      xaxisLabels: model.rLabels,
      yaxisScaleMax: galaxy.mMax,
      xaxisTitle: "Distance (1000 ly)",
      yaxisTitle: "Shell mass (10⁹ M⊙)",
      yaxisTitlePos: 0.3,
      marginLeft: 50,
      marginRight: 10,
      marginTop: 10,
      linewidth: 2,
      colors: ["red", "#ccc"],
      shadow: false,
      dashed: false,
      adjustable: true,
      adjustableOnly: new Array(model.r.length).fill(true),
      outofbounds: true,
      spline: false,
      tickmarksStyle: ["circle", null],
      tickmarksSize: 5
    }
  }).draw();
  densityPlot.on("adjust", function (obj) {
    updateMassModel();
    updateBhSlider();
    updateVelocityPlot();
  });

  if (galaxy.luminousMassModel !== undefined) {
    showLMCheckbox.disabled = false;
    if (showLMCheckbox.checked) {
      densityPlot.show(1);
    } else {
      densityPlot.hide(1);
    }
  } else {
    showLMCheckbox.disabled = true;
    densityPlot.hide(1);
  }

  updateMassModel();
  updateBhSlider();
}

function updateMassModel() {
  var mTotal = model.update(densityPlot.data[0]);
  summaryText.innerHTML =
      "Total mass of galaxy: " +
      formatFriendlyNumber(mTotal * mUnit / MSun) + " solar masses<br/>" +
      "Reduced chi-squared: " + model.calcReducedChiSqrd(galaxy).toFixed(1);
  bhText.innerHTML = formatFriendlyNumber(model.m[0] * mUnit / MSun) +
      " solar masses";
}

function updateVelocityPlot() {
  velocityPlot.data[1] = model.makePlotData();
  RGraph.redraw();
}

function onUpdateBhSlider() {
  var x = bhSlider.value, m = 0;
  if (x > 5) {
    m = Math.pow(10, x / 100) * 100000 * MSun;
  }
  densityPlot.original_data[0][0] = m / mUnit;
  updateMassModel();
  updateVelocityPlot();
}

function updateBhSlider() {
  var x = 0, m = model.m[0] * mUnit;
  if (m > 0) {
    x = 100 * Math.log10(m / (100000 * MSun));
  }
  bhSlider.value = x;
}

function onShowLumMassClicked() {
  if (showLMCheckbox.checked && galaxy.luminousMassModel !== undefined) {
    densityPlot.show(1);
  } else {
    densityPlot.hide(1);
  }
}

function formatFriendlyNumber(x) {
  if (x < 10) {
    return x.toFixed(1);
  }
  if (x < 1000) {
    return x.toFixed();
  }
  if (x < 10000) {
    return (x / 1000).toFixed(1).replace(".", ",") + "00";
  }
  if (x < 1e6) {
    return (x / 1000).toFixed() + ",000";
  }
  const scales = [
    {value: 1e6, name: " million"},
    {value: 1e9, name: " billion"},
    {value: 1e12, name: " trillion"},
  ];
  var i;
  for (i = 0; i < scales.length - 1 && x > scales[i + 1].value; i++);
  x /= scales[i].value;
  return x.toFixed((x < 10) ? 1 : 0) + scales[i].name;
}
