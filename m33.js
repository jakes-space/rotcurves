const M33 = {
  name: "M33",

  // Plot and model settings:
  rMin:  350 * ly,
  rMax:  30 * 1000 * ly,
  rStep: 2 * 1000 * ly,
  vMax:  150 * 1000, // m/s
  mMax:  5 * 1e9 * MSun,

  luminousMassModel: { m0: 1.6, r0: 9.0, mu: 1.0, nu: 1.0 },

  // Rotation curve:
  rUnit: kpc,
  vUnit: 1000, // m/s
  data: [
    [ 0.123 ,   18.1 ,   1.45 ],
    [ 0.245 ,   25.2 ,   0.31 ],
    [ 0.368 ,   31.7 ,   0.23 ],
    [ 0.49  ,   38.7 ,   1.16 ],
    [ 0.613 ,   44.7 ,   1.84 ],
    [ 0.736 ,   49.7 ,   1.84 ],
    [ 0.858 ,   51.9 ,   1.66 ],
    [ 0.981 ,   54.1 ,   0.22 ],
    [ 1.10  ,   56.4 ,   0.81 ],
    [ 1.23  ,   57.9 ,   0.30 ],
    [ 1.35  ,   59.1 ,   0.10 ],
    [ 1.47  ,   62.6 ,   0.64 ],
    [ 1.59  ,   67.3 ,   0.74 ],
    [ 1.72  ,   70.8 ,   0.03 ],
    [ 1.84  ,   73.7 ,   1.23 ],
    [ 1.96  ,   75.9 ,   2.05 ],
    [ 2.08  ,   78.2 ,   1.47 ],
    [ 2.21  ,   79.4 ,   1.44 ],
    [ 2.33  ,   79.8 ,   1.01 ],
    [ 2.45  ,   81.5 ,   0.04 ],
    [ 2.57  ,   83.1 ,   0.58 ],
    [ 2.70  ,   86.2 ,   0.02 ],
    [ 2.82  ,   88.2 ,   0.43 ],
    [ 2.94  ,   88.1 ,   0.67 ],
    [ 3.07  ,   88.7 ,   2.43 ],
    [ 3.19  ,   90.3 ,   3.78 ],
    [ 3.31  ,   92.5 ,   4.00 ],
    [ 3.43  ,   91.8 ,   3.97 ],
    [ 3.56  ,   92.0 ,   3.92 ],
    [ 3.68  ,   91.6 ,   3.81 ],
    [ 3.80  ,   91.4 ,   4.04 ],
    [ 3.92  ,   91.9 ,   4.08 ],
    [ 4.05  ,   93.0 ,   3.17 ],
    [ 4.17  ,   94.6 ,   2.67 ],
    [ 4.29  ,   95.4 ,   3.07 ],
    [ 4.41  ,   95.8 ,   3.82 ],
    [ 4.54  ,   95.3 ,   3.67 ],
    [ 4.66  ,   95.1 ,   3.14 ],
    [ 4.78  ,   96.7 ,   2.45 ],
    [ 4.90  ,   98.1 ,   1.81 ],
    [ 5.03  ,   98.9 ,   0.81 ],
    [ 5.15  ,   99.5 ,   0.06 ],
    [ 5.27  ,   99.3 ,   0.46 ],
    [ 5.39  ,   99.4 ,   0.63 ],
    [ 5.52  ,  101   ,   1.18 ],
    [ 5.64  ,  102   ,   1.60 ],
    [ 5.76  ,  103   ,   2.06 ],
    [ 5.89  ,  104   ,   1.73 ],
    [ 6.01  ,  104   ,   1.42 ],
    [ 6.13  ,  103   ,   1.67 ],
    [ 6.25  ,  103   ,   1.50 ],
    [ 6.38  ,  102   ,   1.11 ],
    [ 6.50  ,  102   ,   1.16 ],
    [ 6.62  ,  103   ,   1.36 ],
    [ 6.74  ,  104   ,   1.37 ],
    [ 6.87  ,  104   ,   1.16 ],
    [ 6.99  ,  104   ,   1.08 ],
    [ 7.11  ,  104   ,   1.06 ],
    [ 7.23  ,  104   ,   1.16 ],
    [ 7.36  ,  105   ,   1.38 ],
    [ 7.48  ,  106   ,   2.33 ],
    [ 7.60  ,  106   ,   2.08 ],
    [ 7.72  ,  106   ,   2.27 ],
    [ 7.85  ,  107   ,   2.85 ],
    [ 7.97  ,  108   ,   2.96 ],
    [ 8.09  ,  107   ,   2.03 ],
    [ 8.21  ,  107   ,   1.13 ],
    [ 8.34  ,  107   ,   0.83 ],
    [ 8.46  ,  111   ,   0.31 ],
    [ 8.58  ,  117   ,   3.14 ],
    [ 8.71  ,  116   ,   3.00 ],
    [ 8.83  ,  118   ,   3.34 ],
    [ 8.95  ,  115   ,   0.12 ],
  ]
}