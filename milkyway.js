const MilkyWay = {
  name: "Milky Way",

  // Plot and model settings:
  rMin:  500 * ly,
  rMax:  55 * 1000 * ly,
  rStep: 5 * 1000 * ly,
  vMax:  400 * 1000, // m/s
  mMax:  100 * 1e9 * MSun,

  specialPoints: [
    { r: 28, v: 240, ptRad: 5, color: "yellow" }  // Sun
  ],

//  luminousMassModel: { m0: 55, r0: 6.5, mu: 1.0, nu: 0.75 },

  // Rotation curve from Bhattacharjee, P., Chaudhury, S., and Susmita, K.
  // (ApJ 2014; https://iopscience.iop.org/article/10.1088/0004-637X/785/1/63)
  // We use the data based on an 8.3 kpc Sun-GC distance.
  rUnit: kpc,
  vUnit: 1000, // m/s
  data: [
    [   0.20 ,  233.0  ,  13.32 ],
    [   0.38 ,  268.92 ,   4.67 ],
    [   0.66 ,  250.75 ,  11.35 ],
    [   1.61 ,  217.83 ,   5.81 ],
    [   2.57 ,  219.58 ,   1.48 ],
    [   3.59 ,  223.11 ,   2.43 ],
    [   4.51 ,  247.88 ,   2.99 ],
    [   5.53 ,  253.14 ,   1.69 ],
    [   6.50 ,  270.95 ,   2.19 ],
    [   7.56 ,  267.80 ,   0.96 ],
    [   8.34 ,  270.52 ,   0.66 ],
    [   9.45 ,  235.58 ,   8.44 ],
    [  10.50 ,  249.72 ,  13.44 ],
    [  11.44 ,  261.96 ,  11.71 ],
    [  12.51 ,  284.30 ,  17.50 ],
    [  13.53 ,  271.54 ,  15.57 ],
    [  14.59 ,  251.43 ,  25.60 ],
    [  16.05 ,  320.70 ,  25.27 ],  // Set rMax = 55 kly here

//    [  18.64 ,  286.46 , 101.18 ],
//    [  26.30 ,  189.64 ,   6.74 ],
//    [  28.26 ,  237.99 ,  11.54 ],
//    [  29.51 ,  209.82 ,   9.16 ],
//    [  32.04 ,  179.14 ,   6.65 ],
//    [  33.99 ,  170.37 ,   6.93 ],
//    [  36.49 ,  175.92 ,   6.62 ],
//    [  38.41 ,  191.57 ,  11.73 ],
//    [  40.42 ,  197.59 ,  14.12 ],
//    [  42.40 ,  192.79 ,   5.92 ],
//    [  44.49 ,  213.22 ,  17.17 ],
//    [  45.99 ,  179.39 ,  11.23 ],
//    [  48.06 ,  213.03 ,  24.72 ],
//    [  49.49 ,  178.57 ,  17.63 ],
//    [  51.39 ,  183.31 ,  23.58 ],
//    [  53.89 ,  157.89 ,  19.57 ],
//    [  56.89 ,  191.76 ,  24.35 ],
//    [  57.98 ,  210.72 ,  29.81 ],
//    [  60.92 ,  168.02 ,  25.67 ],  // Set rMax = 200 kly here

//    [  64.73 ,  206.47 ,  36.27 ],
//    [  69.31 ,  203.62 ,  40.89 ],
//    [  72.96 ,  190.53 ,  40.98 ],
//    [  76.95 ,  222.72 ,  74.37 ],
//    [  81.13 ,  186.29 ,  66.53 ],
//    [  84.90 ,  122.25 ,  36.46 ],
//    [  89.35 ,  143.95 ,  29.49 ],
//    [  92.44 ,  154.66 ,  67.23 ],
//    [  97.41 ,  184.0  ,  72.86 ],
//    [ 100.72 ,  108.68 ,  40.99 ],
//    [ 106.77 ,  137.15 ,  53.17 ],
//    [ 119.98 ,  150.18 ,  25.46 ],
//    [ 189.49 ,  125.01 ,  37.32 ],
  ]
}
