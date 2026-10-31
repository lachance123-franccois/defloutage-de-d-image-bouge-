# 18 – Défloutage de bougé : estimation du noyau, Wiener, Richardson-Lucy, variation totale

Restaurer une photo floue à cause d'un mouvement de l'appareil pendant la pose. Le flou est modélisé par une convolution : y = k ∗ x + n, où k est le **noyau de flou** (la trajectoire de l'appareil). Le projet traite les deux questions :

1. **Déconvolution à noyau connu** : Wiener, Richardson-Lucy, régularisation par variation totale (TV).
2. **Estimation du noyau à partir de la seule photo floue** : par le cepstre pour un bougé rectiligne, et par une méthode aveugle multi-échelle pour un tremblement quelconque.

Tout est écrit **en MATLAB, sans toolbox**.

## Lancer

```matlab
cd 18-defloutage-bouge
main_nonblind              % noyau connu : 3 méthodes, 4 images, 4 noyaux (~5 min)
main_kernel_estimation     % cepstre + estimation aveugle (~5 min)
deblur_photo               % votre photo : data/photo.jpg (sinon démonstration)
```

## Fichiers

| Fichier | Rôle |
|---|---|
| `main_nonblind.m` | benchmark à noyau connu, semi-convergence de RL, choix de λ |
| `main_kernel_estimation.m` | estimation du noyau (cepstre, aveugle) puis déconvolution |
| `deblur_photo.m` | chaîne complète sur une photo couleur |
| `motion_kernel.m` | bougé rectiligne anti-aliasé (équivalent de `fspecial('motion')`) |
| `shake_kernel.m`, `splat_points.m` | tremblement aléatoire (trajectoire lisse à vitesse variable) |
| `blur_image.m` | flou réaliste, sans repliement circulaire |
| `deconv_wiener.m`, `edge_taper.m`, `psf_to_otf.m` | filtre de Wiener par FFT, traitement des bords |
| `deconv_rl.m` | Richardson-Lucy avec modèle de bords correct |
| `deconv_tv.m` | variation totale, algorithme primal-dual de Chambolle-Pock |
| `conv_fft.m` | convolutions 'valid' et 'full' par FFT |
| `estimate_motion_cepstrum.m` | longueur et angle d'un bougé rectiligne |
| `blind_kernel_estimate.m` | estimation aveugle (prior L0, grossier → fin) |
| `psnr_aligned.m`, `kernel_similarity.m` | évaluation des méthodes aveugles |
| `data/kernels.mat` | les 4 noyaux du benchmark (K1 à K4) |

## Méthodes

**Modèle de flou sans repliement.** Une convolution par FFT suppose l'image périodique : le bord droit « déborde » sur le bord gauche. Une vraie photo ne fonctionne pas ainsi, la scène continue hors du cadre. Le flou est donc simulé par une extension des bords suivie d'une convolution `'valid'`, et RL et TV utilisent le même modèle : l'inconnue x est plus grande que y d'une demi-taille de noyau de chaque côté (Bertero & Boccacci, 2005). Il n'y a ainsi aucune oscillation de bord à corriger après coup.

**Wiener.** X(ω) = K̄(ω) Y(ω) / (|K(ω)|² + NSR). Là où le spectre du flou s'annule, le terme NSR empêche d'amplifier le bruit. Pour la FFT, l'image est étendue par symétrie puis atténuée vers les bords (`edge_taper`).

**Richardson-Lucy.** Maximum de vraisemblance pour un bruit de Poisson, avec des itérations multiplicatives qui gardent x ≥ 0 : x ← x · Aᵀ(y / Ax) / Aᵀ1. Le nombre d'itérations sert de régularisation. Au-delà d'environ 40 itérations le bruit est amplifié et le PSNR redescend (*semi-convergence*, voir `results/parameters.png`).

**Variation totale.** min ½‖Ax − y‖² + λ TV(x), avec 0 ≤ x ≤ 1, résolu par l'algorithme primal-dual de Chambolle & Pock (2011). La TV favorise les zones lisses séparées de contours nets, ce qui supprime à la fois le bruit et les oscillations. λ ≈ 0,1 × σ du bruit.

**Cepstre (bougé rectiligne).** Le spectre d'un bougé de longueur L est un sinus cardinal, qui s'annule périodiquement : la photo floue montre des bandes sombres perpendiculaires au mouvement. Le cepstre ifft(log|Y|) transforme cette périodicité en deux pics négatifs, situés à ±L pixels dans la direction du mouvement. On lit donc L et l'angle sur la position du minimum.

**Estimation aveugle (noyau quelconque)**, d'après Xu, Zheng & Jia (2013) et Pan et al. (2014). On alterne :
- image latente : min ‖k∗x − y‖² + λ‖∇x‖₀. Le prior L0 ne garde que les contours francs, les plus informatifs sur le flou ;
- noyau : min ‖k∗∇x − ∇y‖² + γ‖k‖², en forme fermée par FFT, puis positivité, seuillage, suppression des petites composantes isolées et recentrage.

Le tout se fait du grossier au fin (pyramide de facteur √2), ce qui évite les minima locaux.

## Résultats de référence

Valeurs obtenues pendant la mise au point avec une implémentation de référence des mêmes algorithmes et paramètres. Les chiffres MATLAB doivent être très proches, la réalisation du bruit étant différente. Bruit gaussien σ = 0,01, images 256×256.

**Noyau connu**, PSNR moyen sur les 4 images (dB) :

| Noyau | Flou | Wiener | Richardson-Lucy (40 it.) | TV (λ = 10⁻³) |
|---|---|---|---|---|
| K1 linéaire 15 px, 30° | 21,40 | 24,54 | 25,51 | **27,10** |
| K2 linéaire 25 px, 120° | 19,94 | 22,88 | 24,01 | **25,27** |
| K3 tremblement (boucle) | 18,71 | 24,66 | 25,98 | **28,65** |
| K4 tremblement | 18,97 | 26,73 | 27,52 | **30,90** |

SSIM moyen, de l'image floue à la TV : de 0,49 à 0,90 pour K4, de 0,60 à 0,81 pour K1.

**Cepstre** (bougé réel K1 = 15 px à 30°, K2 = 25 px à 120°) : estimations de 15,3 px / 31,6° pour K1 et de 23,7 à 24,2 px / 117,6° à 119,7° pour K2 sur cameraman, astronaut et coffee. La texture périodique de « brick » trompe la méthode (14,3 px / 24,8° pour K1), car ses propres répétitions créent des pics dans le cepstre. Avec le noyau estimé, la TV perd de 0,04 à 0,8 dB par rapport au vrai noyau (hors brick).

**Estimation aveugle** : similarité entre noyau estimé et vrai noyau (1 = identiques), puis PSNR de la TV avec le noyau estimé, contre la TV avec le vrai noyau :

| Noyau | Similarité (camera / astronaut / coffee) | TV noyau estimé | TV vrai noyau |
|---|---|---|---|
| K1 | 0,92 / 0,91 / 0,93 | 25,1 / 24,0 / 26,2 dB | 26,0 / 25,5 / 27,6 dB |
| K2 | 0,90 / 0,90 / 0,90 | 22,7 / 22,1 / 24,2 dB | 24,2 / 23,9 / 25,1 dB |
| K3 | 0,85 / 0,79 / 0,71 | 26,0 / 24,6 / 24,8 dB | 27,5 / 27,3 / 29,0 dB |
| K4 | 0,91 / 0,91 / 0,88 | 29,1 / 28,1 / 28,2 dB | 30,9 / 29,2 / 29,8 dB |

L'estimation aveugle récupère l'essentiel du gain : de 18 dB (flou) à 28-29 dB pour K4, soit 1 à 2 dB sous la déconvolution avec le vrai noyau. Elle échoue sur « brick » (similarité 0,3 à 0,6) : une texture répétitive sans contour dominant ne contraint pas le noyau.

## Limites et pistes

- Flou supposé uniforme sur toute l'image (pas de rotation de l'appareil, pas d'objet en mouvement).
- Les PSNR des méthodes aveugles sont mesurés au meilleur recalage (±8 px), le noyau n'étant défini qu'à une translation près.
- Pistes : prior sur les intensités sombres (Pan 2016), noyaux variables dans l'image, réseaux de déconvolution déroulés (dans la continuité du projet *deconvolution-variational-to-unrolled*).

## Références

- N. Wiener, *Extrapolation, Interpolation, and Smoothing of Stationary Time Series*, 1949.
- W. H. Richardson (1972), L. B. Lucy (1974) ; M. Bertero, P. Boccacci, *A simple method for the reduction of boundary effects in the Richardson-Lucy approach*, A&A 2005.
- L. Rudin, S. Osher, E. Fatemi, *Nonlinear total variation based noise removal*, Physica D 1992.
- A. Chambolle, T. Pock, *A first-order primal-dual algorithm for convex problems*, JMIV 2011.
- M. Cannon, *Blind deconvolution of spatially invariant image blurs with phase*, IEEE TASSP 1976 ; R. Rom, *On the cepstrum of two-dimensional functions*, IEEE TIT 1975.
- L. Xu, S. Zheng, J. Jia, *Unnatural L0 sparse representation for natural image deblurring*, CVPR 2013.
- J. Pan, Z. Hu, Z. Su, M.-H. Yang, *Deblurring text images via L0-regularized intensity and gradient prior*, CVPR 2014.
- A. Levin, Y. Weiss, F. Durand, W. Freeman, *Understanding and evaluating blind deconvolution algorithms*, CVPR 2009.
- G. Boracchi, A. Foi, *Modeling the performance of image restoration from motion blur*, IEEE TIP 2012.
