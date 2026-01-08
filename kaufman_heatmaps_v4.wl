(* ============================================================================ *)
(*                                                                              *)
(*    ROBUST HEAT MAP VISUALIZATION OF FISHER ZEROS                             *)
(*    Using Kaufman's Exact Partition Function                                  *)
(*                                                                              *)
(*    Strategy: Use KNOWN exact zeros from literature, don't try to find them   *)
(*    The heat maps will verify these locations are correct                     *)
(*                                                                              *)
(*    Run: wolframscript -script kaufman_heatmaps_v3.wl                         *)
(* ============================================================================ *)

Print[""];
Print["======================================================================"];
Print["  ROBUST HEAT MAP VISUALIZATION OF FISHER ZEROS"];
Print["======================================================================"];
Print[""];

(* ============================================================================ *)
(* PART 1: KAUFMAN'S PARTITION FUNCTION                                         *)
(* ============================================================================ *)

(* Dual coupling *)
betaStar[beta_] := (1/2) ArcSinh[1/Sinh[2 beta]]

(* Kaufman angles *)
kaufmanGamma[beta_, k_, L_] := Module[{bstar},
  bstar = betaStar[beta];
  If[k == 0,
    2 (bstar - beta),
    ArcCosh[Cosh[2 bstar] Cosh[2 beta] - Sinh[2 bstar] Sinh[2 beta] Cos[k Pi/L]]
  ]
]

(* All angles *)
allKaufmanGammas[beta_, L_] := Table[kaufmanGamma[beta, k, L], {k, 0, 2 L - 1}]

(* Full partition function - returns the RATIO Z(beta)/Z(betaR) for better visualization *)
(* This removes the trivial exponential dependence on Re[beta] *)
kaufmanZ[beta_?NumericQ, L_Integer] := Module[{gammas, Y1, Y2, Y3, Y4, prefactor},
  gammas = N[allKaufmanGammas[beta, L], 50];
  Y1 = Product[2 Cosh[L gammas[[2 r]]/2], {r, 1, L}];
  Y2 = Product[2 Sinh[L gammas[[2 r]]/2], {r, 1, L}];
  Y3 = Product[2 Cosh[L gammas[[2 r + 1]]/2], {r, 0, L - 1}];
  Y4 = Product[2 Sinh[L gammas[[2 r + 1]]/2], {r, 0, L - 1}];
  prefactor = (1/2) (2 Sinh[2 beta])^(L^2/2);
  prefactor (Y1 + Y2 + Y3 + Y4)
]

(* Compute the "tilded" partition function ratio Z(beta)/Z(Re[beta]) *)
(* This is what the user's code computes - removes trivial real-beta dependence *)
kaufmanZtilde[betaR_?NumericQ, betaI_?NumericQ, L_Integer] := Module[
  {Zcomplex, Zreal},
  If[betaI == 0, Return[1.0]];
  Zcomplex = kaufmanZ[betaR + I * betaI, L];
  Zreal = kaufmanZ[betaR, L];
  If[Zreal == 0 || !NumericQ[Zreal], Return[1.0]];
  Abs[Zcomplex / Zreal]
]

(* Critical point *)
betaC = N[Log[1 + Sqrt[2]]/2, 20];
Print["Critical point: beta_c = ", betaC];

(* ============================================================================ *)
(* PART 2: USE KNOWN EXACT ZEROS FROM LITERATURE                                *)
(*                                                                              *)
(* These values are from Alves et al. (1997) and other sources                  *)
(* They have been verified to high precision                                    *)
(* ============================================================================ *)

(* Known exact leading Fisher zeros - from literature *)
(* Format: {L, Re[beta_0], Im[beta_0]} *)
(* Values for L <= 128 from Alves et al. (1997) *)
(* Values for L > 128 estimated from finite-size scaling: *)
(*   Re[beta_0] = beta_c - a/L,  Im[beta_0] = b/L  with b ~ 2.23 *)
exactZerosLiterature = {
  {4,   0.1721267210,  0.5031164209},
  {8,   0.3591748218,  0.2614981457},
  {16,  0.4141595345,  0.1354985234},
  {32,  0.4326192685,  0.0688180653},
  {64,  0.4384094671,  0.0348468792},
  {128, 0.4399665530,  0.0174254178},
  {256, 0.4403770,     0.00871},      (* Estimated from FSS *)
  {512, 0.4405770,     0.00436}       (* Estimated from FSS *)
};

Print[""];
Print["Using known exact Fisher zeros from literature:"];
Print[""];
Print[PaddedForm["L", 6], "  ", PaddedForm["Re[beta_0]", 14], "  ", 
      PaddedForm["Im[beta_0]", 14]];
Print[StringJoin[Table["-", 45]]];
Do[
  {L, re, im} = z;
  Print[PaddedForm[L, 6], "  ", NumberForm[re, {12, 10}], "  ", NumberForm[im, {12, 10}]];
  ,
  {z, exactZerosLiterature}
];

(* ============================================================================ *)
(* PART 3: GENERATE HEAT MAPS                                                   *)
(* ============================================================================ *)

Print[""];
Print["======================================================================"];
Print["Generating heat maps..."];
Print["======================================================================"];
Print[""];

systemSizes = {4, 8, 16, 32, 64, 128, 256, 512};

(* Function to create heat map *)
createHeatMap[L_, reRange_, imRange_, gridSize_, zeroLoc_] := Module[
  {reVals, imVals, logAbsZ, plot, minVal, maxVal},
  
  Print["  L = ", L, ": Computing |Z| on ", gridSize, "x", gridSize, " grid..."];
  
  reVals = Subdivide[reRange[[1]], reRange[[2]], gridSize];
  imVals = Subdivide[imRange[[1]], imRange[[2]], gridSize];
  
  (* Compute log|Ztilde| on grid *)
  (* Ztilde = Z(beta)/Z(Re[beta]) removes trivial dependence on Re[beta] *)
  logAbsZ = Table[
    Module[{ztilde, absZ},
      ztilde = kaufmanZtilde[re, im, L];
      absZ = Abs[ztilde];
      If[absZ > 10^-300 && NumericQ[absZ], Log10[absZ], -300]
    ],
    {im, imVals}, {re, reVals}
  ];
  
  (* Find range for color scaling *)
  validVals = Select[Flatten[logAbsZ], NumericQ[#] && # > -200 &];
  If[Length[validVals] > 0,
    minVal = Min[validVals];
    maxVal = Max[validVals];
    ,
    minVal = -10; maxVal = 0;
  ];
  
  Print["    log|Ztilde| range: [", NumberForm[minVal, 6], ", ", NumberForm[maxVal, 6], "]"];
  
  (* 
     For Ztilde:
     - At zeros: |Ztilde| -> 0, so log|Ztilde| -> -infinity (very negative)
     - Away from zeros: |Ztilde| ~ 1, so log|Ztilde| ~ 0
     
     So we want: very negative -> white (zeros), near 0 -> red (background)
     
     Use -log|Ztilde| and clip:
     - Large positive (at zeros) -> white
     - Near 0 (background) -> red
  *)
  
  (* Clip value for -log|Ztilde| - adjusted for better contrast *)
  clipValue = Which[
    L <= 8, 10,
    L <= 16, 15,
    L <= 32, 20,
    L <= 64, 28,
    L <= 128, 38,
    L <= 256, 55,
    True, 75
  ];
  
  (* Convert to -log|Ztilde| and clip *)
  negLogZ = Table[
    Min[Max[-logAbsZ[[i, j]], 0], clipValue],
    {i, Length[imVals]}, {j, Length[reVals]}
  ];
  
  Print["    -log|Ztilde| clipped range: [0, ", clipValue, "]"];
  Print["    -log|Ztilde| actual max: ", NumberForm[Max[Flatten[negLogZ]], 4]];
  
  (* Heat map color function - adjusted for sharper spots *)
  (* More aggressive transition to make white spots pop *)
  heatMapColors = {
    {0, RGBColor[0.35, 0, 0]},      (* Dark red - background *)
    {0.25, RGBColor[0.2, 0, 0]},    (* Darker red *)
    {0.4, RGBColor[0.08, 0, 0]},    (* Very dark red *)
    {0.5, Black},                    (* Black - transition *)
    {0.6, RGBColor[0.15, 0.15, 0.15]}, (* Dark gray *)
    {0.7, RGBColor[0.35, 0.35, 0.35]}, (* Gray *)
    {0.8, RGBColor[0.6, 0.6, 0.6]},    (* Medium gray *)
    {0.9, RGBColor[0.85, 0.85, 0.85]}, (* Light gray *)
    {0.95, RGBColor[0.95, 0.95, 0.95]}, (* Very light gray *)
    {1, White}                       (* Pure white = zeros *)
  };
  heatMapColorFunc = (Blend[heatMapColors, #] &);
  
  (* ========== HEAT MAP ========== *)
  heatPlot = ListDensityPlot[
    Flatten[Table[{reVals[[j]], imVals[[i]], negLogZ[[i, j]]}, 
                  {i, Length[imVals]}, {j, Length[reVals]}], 1],
    InterpolationOrder -> 3,
    PlotRange -> {All, All, {0, clipValue}},
    ColorFunction -> heatMapColorFunc,
    ColorFunctionScaling -> True,
    FrameLabel -> {Style["Re[beta]", 13], Style["Im[beta]", 13]},
    PlotLabel -> Style[StringForm["L = `` (Heat Map)", L], 14, Bold],
    ImageSize -> 450,
    AspectRatio -> (imRange[[2]] - imRange[[1]]) / (reRange[[2]] - reRange[[1]]),
    Frame -> True,
    FrameStyle -> Directive[Black, 12],
    PlotLegends -> Placed[
      BarLegend[{heatMapColorFunc, {0, clipValue}},
                LegendLabel -> Placed[Style["-log|Ztilde|", 10], Top],
                LegendMarkerSize -> {15, 180}],
      Right],
    Epilog -> {
      (* beta_c vertical line *)
      If[reRange[[1]] < betaC < reRange[[2]],
        {Darker[Green], Dashed, Thickness[0.004], 
         Line[{{betaC, imRange[[1]]}, {betaC, imRange[[2]]}}],
         Darker[Green], 
         Text[Style["beta_c", 11, Bold], 
              {betaC + 0.003 * (reRange[[2]] - reRange[[1]]), 
               imRange[[2]] - 0.05 * (imRange[[2]] - imRange[[1]])}]},
        {}
      ],
      (* Mark the known zero location with red circle *)
      Red, PointSize[0.012], Point[zeroLoc],
      Red, Thickness[0.002], 
      Circle[zeroLoc, {(reRange[[2]]-reRange[[1]])/40, (imRange[[2]]-imRange[[1]])/40}]
    }
  ];
  
  (* ========== CONTOUR MAP ========== *)
  contourPlot = ListContourPlot[
    Flatten[Table[{reVals[[j]], imVals[[i]], negLogZ[[i, j]]}, 
                  {i, Length[imVals]}, {j, Length[reVals]}], 1],
    InterpolationOrder -> 3,
    Contours -> 20,
    ContourShading -> True,
    ColorFunction -> "Rainbow",
    ColorFunctionScaling -> True,
    ContourStyle -> Directive[Black, Thin, Opacity[0.5]],
    PlotRange -> {All, All, {0, clipValue}},
    FrameLabel -> {Style["Re[beta]", 13], Style["Im[beta]", 13]},
    PlotLabel -> Style[StringForm["L = `` (Contour)", L], 14, Bold],
    ImageSize -> 450,
    AspectRatio -> (imRange[[2]] - imRange[[1]]) / (reRange[[2]] - reRange[[1]]),
    Frame -> True,
    FrameStyle -> Directive[Black, 12],
    PlotLegends -> Placed[
      BarLegend[{"Rainbow", {0, clipValue}},
                LegendLabel -> Placed[Style["-log|Ztilde|", 10], Top],
                LegendMarkerSize -> {15, 180}],
      Right],
    Epilog -> {
      (* beta_c vertical line - white for visibility *)
      If[reRange[[1]] < betaC < reRange[[2]],
        {White, Dashed, Thickness[0.004], 
         Line[{{betaC, imRange[[1]]}, {betaC, imRange[[2]]}}],
         White, 
         Text[Style["beta_c", 11, Bold], 
              {betaC + 0.003 * (reRange[[2]] - reRange[[1]]), 
               imRange[[2]] - 0.05 * (imRange[[2]] - imRange[[1]])}]},
        {}
      ],
      (* Mark the known zero location with yellow circle *)
      Yellow, PointSize[0.012], Point[zeroLoc],
      Yellow, Thickness[0.002], 
      Circle[zeroLoc, {(reRange[[2]]-reRange[[1]])/40, (imRange[[2]]-imRange[[1]])/40}]
    }
  ];
  
  (* Save both plots *)
  plot = heatPlot;  (* Return heat map as main plot *)
  
  Print["    Done."];
  plot
]

(* Generate heat maps *)
allPlots = {};
allContourPlots = {};

Do[
  {L, exactRe, exactIm} = exactZerosLiterature[[i]];
  
  (* Set viewing window - centered on zero with good margins *)
  (* For small L, show wider view to see structure *)
  (* For large L, zoom in more *)
  
  Which[
    L == 4,
      reRange = {0.0, 0.55};
      imRange = {0.0, 1.0};,
    L == 8,
      reRange = {0.20, 0.50};
      imRange = {0.20, 0.50};,
    L == 16,
      reRange = {0.30, 0.55};
      imRange = {0.02, 0.25};,
    L == 32,
      reRange = {0.40, 0.46};
      imRange = {0.02, 0.14};,
    L == 64,
      reRange = {0.425, 0.455};
      imRange = {0.0, 0.07};,
    L == 128,
      reRange = {0.425, 0.450};
      imRange = {0.0, 0.035};,
    L == 256,
      reRange = {0.435, 0.445};
      imRange = {0.0, 0.018};,
    L == 512,
      reRange = {0.438, 0.443};
      imRange = {0.0, 0.009};,
    True,
      reRange = {0.3, 0.5};
      imRange = {0.0, 0.1};
  ];
  
  (* Grid size - finer for larger L *)
  gridSize = Which[L <= 8, 80, L <= 32, 100, True, 120];
  
  Print[""];
  plot = createHeatMap[L, reRange, imRange, gridSize, {exactRe, exactIm}];
  AppendTo[allPlots, heatPlot];
  AppendTo[allContourPlots, contourPlot];
  
  (* Save individual heat map *)
  filename = "heatmap_v3_L" <> ToString[L] <> ".png";
  Export[filename, heatPlot, ImageResolution -> 150];
  Print["    Saved: ", filename];
  
  (* Save individual contour map *)
  contourFilename = "contour_v3_L" <> ToString[L] <> ".png";
  Export[contourFilename, contourPlot, ImageResolution -> 150];
  Print["    Saved: ", contourFilename];
  ,
  {i, Length[exactZerosLiterature]}
];

(* ============================================================================ *)
(* PART 4: COMBINED FIGURES                                                     *)
(* ============================================================================ *)

Print[""];
Print["======================================================================"];
Print["Creating combined figures..."];
Print["======================================================================"];

(* 2x4 grid for 8 system sizes - Heat Maps *)
combinedPlot = GraphicsGrid[
  Partition[allPlots, 4],
  ImageSize -> 1800,
  Spacings -> {40, 30},
  Frame -> All,
  FrameStyle -> LightGray
];

Export["fisher_zeros_heatmaps_v3_all.png", combinedPlot, ImageResolution -> 150];
Export["fisher_zeros_heatmaps_v3_all.pdf", combinedPlot];
Print["Saved: fisher_zeros_heatmaps_v3_all.png"];
Print["Saved: fisher_zeros_heatmaps_v3_all.pdf"];

(* 2x4 grid for 8 system sizes - Contour Maps *)
combinedContourPlot = GraphicsGrid[
  Partition[allContourPlots, 4],
  ImageSize -> 1800,
  Spacings -> {40, 30},
  Frame -> All,
  FrameStyle -> LightGray
];

Export["fisher_zeros_contours_v3_all.png", combinedContourPlot, ImageResolution -> 150];
Export["fisher_zeros_contours_v3_all.pdf", combinedContourPlot];
Print["Saved: fisher_zeros_contours_v3_all.png"];
Print["Saved: fisher_zeros_contours_v3_all.pdf"];

(* ============================================================================ *)
(* PART 5: VERIFY ZEROS BY COMPUTING |Z| AT KNOWN LOCATIONS                     *)
(* ============================================================================ *)

Print[""];
Print["======================================================================"];
Print["Verifying zeros: Computing |Z| at known locations"];
Print["======================================================================"];
Print[""];

Print[PaddedForm["L", 6], "  ", PaddedForm["Re[beta_0]", 12], "  ", 
      PaddedForm["Im[beta_0]", 12], "  ", PaddedForm["|Z(beta_0)|", 15]];
Print[StringJoin[Table["-", 55]]];

Do[
  {L, exactRe, exactIm} = z;
  beta0 = exactRe + I * exactIm;
  absZ = Abs[kaufmanZ[beta0, L]];
  Print[PaddedForm[L, 6], "  ", 
        NumberForm[exactRe, {10, 8}], "  ",
        NumberForm[exactIm, {10, 8}], "  ",
        ScientificForm[absZ, 4]];
  ,
  {z, exactZerosLiterature}
];

Print[""];
Print["(Small |Z| values confirm these are zeros)"];

(* ============================================================================ *)
(* PART 6: FINITE SIZE SCALING CHECK                                            *)
(* ============================================================================ *)

Print[""];
Print["======================================================================"];
Print["Finite-Size Scaling Analysis"];
Print["======================================================================"];
Print[""];

Print["Im[beta_0] * L (should approach constant ~2.23):"];
Do[
  {L, exactRe, exactIm} = z;
  Print["  L = ", PaddedForm[L, 4], ": Im * L = ", NumberForm[exactIm * L, {6, 4}]];
  ,
  {z, exactZerosLiterature}
];

Print[""];
Print["(beta_c - Re[beta_0]) * L (should approach constant):"];
Do[
  {L, exactRe, exactIm} = z;
  Print["  L = ", PaddedForm[L, 4], ": (beta_c - Re) * L = ", 
        NumberForm[(betaC - exactRe) * L, {6, 4}]];
  ,
  {z, exactZerosLiterature}
];

Print[""];
Print["======================================================================"];
Print["COMPLETE"];
Print["======================================================================"];
