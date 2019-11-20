const UGC128 = {
  name: "UGC128",

  // Plot and model settings:
  rMin:  1000 * ly,
  rMax:  150 * 1000 * ly,
  rStep: 25 * 1000 * ly,
  vMax:  150 * 1000, // m/s
  mMax:  50 * 1e9 * MSun,

  // This should provide only about 5% of the galaxy's mass.
  luminousMassModel: { m0: 14.5, r0: 6.5, mu: 1.0, nu: 0.66 },

  // Rotation curve:
  rUnit: kpc,
  vUnit: 1000, // m/s
  data: [
    [  1.6 ,   29.3 ,  8 ],
    [  3.1 ,   51.3 ,  8 ],
    [  4.6 ,   65.9 ,  8 ],
    [  6.3 ,   78.1 ,  8 ],
    [  9.3 ,   94.0 ,  8 ],
    [ 12.4 ,  107.4 ,  8 ],
    [ 15.5 ,  117.2 ,  8 ],
    [ 21.7 ,  125.7 ,  8 ],
    [ 27.9 ,  129.4 ,  8 ],
    [ 34.1 ,  130.6 ,  8 ],
    [ 40.4 ,  130.6 ,  8 ],
    [ 45.3 ,  130.6 ,  8 ],
  ]
}
