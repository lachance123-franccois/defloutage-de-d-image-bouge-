# Défloutage de bougé

Restaurer une photo floue à cause d'un mouvement de l'appareil pendant la pose : estimation du flou, filtre de Wiener, Richardson-Lucy, variation totale.

## Le problème

Quand l'appareil bouge pendant la pose, chaque point de la scène laisse une traînée au lieu d'un point net. La photo floue se modélise par

$$y = k * x + n$$

- $x$ : l'image nette, inconnue
- $k$ : le **noyau de flou**, une petite image qui dessine la trajectoire de l'appareil ($k \ge 0$, $\sum k = 1$)
- $n$ : un bruit gaussien d'écart-type $\sigma = 0{,}01$
- $y$ : la photo observée
- $*$ : la convolution, $(k*x)(i,j)=\sum_{u,v} k(u,v)\,x(i-u,\,j-v)$

> **En clair :** chaque pixel de la photo floue est une moyenne de ses voisins, pondérée par la forme de la trajectoire de l'appareil.

Deux questions :

1. Si $k$ est connu, comment retrouver $x$ ? (déconvolution non aveugle)
2. Si $k$ est inconnu, comment le deviner d'abord ? (déconvolution aveugle)

---

## 1. Simuler un flou réaliste

### Les noyaux

**Bougé rectiligne.** Un segment de longueur $L$ et d'angle $\theta$ :

$$p(t) = t\,(\cos\theta,\ -\sin\theta), \qquad t \in \left[-\tfrac{L}{2}, \tfrac{L}{2}\right]$$

Le segment est échantillonné finement, et chaque point (à coordonnées non entières) est réparti sur les 4 pixels voisins avec des poids bilinéaires. Cela évite l'effet d'escalier.

**Tremblement aléatoire.** Une courbe lisse, somme de 3 sinusoïdes lentes :

$$p_x(t)=\sum_{j=1}^{3} a_j \sin(2\pi f_j t+\varphi_j), \qquad p_y(t)=\sum_{j=1}^{3} b_j \sin(2\pi f'_j t+\varphi'_j)$$

avec $f_j \in [0{,}2 ;\ 1{,}2]$. Elle est parcourue à vitesse variable : là où l'appareil ralentit, il reste plus longtemps, donc le noyau y est plus lumineux.

### Les bords

Une convolution par FFT suppose l'image périodique, comme un carrelage : le bord droit « déborde » sur le bord gauche. Une vraie photo ne fonctionne pas ainsi, la scène continue hors du cadre.

On prolonge donc l'image en répétant ses pixels de bord, puis on convolue en mode « valide » (sans rien inventer au-delà). L'image nette inconnue est alors plus grande que la photo, d'une demi-taille de noyau de chaque côté. Cela supprime les oscillations parasites aux bords.

---

## 2. Déconvolution quand le noyau est connu

### Wiener

Dans Fourier, la convolution devient une multiplication : $Y = K \cdot X$. Diviser par $K$ est impossible là où $|K| \approx 0$ : le bruit y serait multiplié par un nombre énorme. Wiener ajoute une constante :

$$\hat X(\omega)=\frac{\overline{K(\omega)}\,Y(\omega)}{|K(\omega)|^2+\mathrm{NSR}}$$

Pour limiter les oscillations dues à la FFT, l'image est étendue par symétrie, puis estompée près des bords en la mélangeant avec sa version floutée circulairement.

> **En clair :** on défait le flou fréquence par fréquence, mais on s'interdit d'amplifier ce qui est trop faible pour être fiable. NSR règle ce niveau de prudence.

### Richardson-Lucy

Pour un bruit de Poisson (comptage de photons), le maximum de vraisemblance conduit à des itérations multiplicatives :

$$x^{(t+1)} = x^{(t)} \cdot \frac{A^{\top}\!\left(y \,/\, A x^{(t)}\right)}{A^{\top}\mathbf{1}}$$

où $A$ est l'opérateur de flou (convolution valide) et $A^{\top}$ son adjoint (corrélation).

> **En clair :** on floute l'image courante, on la compare à la photo réelle par un rapport, puis on corrige l'image avec ce rapport. Comme on multiplie, l'image reste positive.

Le nombre d'itérations sert de régularisation : au-delà d'une quarantaine, l'algorithme se met à reproduire le bruit et la qualité redescend. C'est la **semi-convergence**.

### Variation totale (TV)

On cherche l'image qui colle aux données **et** qui reste simple :

$$\min_{0 \le x \le 1}\ \tfrac12 \lVert Ax-y\rVert^2 + \lambda\,\mathrm{TV}(x), \qquad \mathrm{TV}(x)=\sum_{i,j}\lVert \nabla x(i,j)\rVert_2$$

> **En clair :** une photo naturelle est faite de zones presque uniformes séparées par des contours nets. Le terme TV pénalise le bruit, qui crée des variations partout, mais tolère quelques contours francs. Le paramètre $\lambda$ dose ce compromis : on prend $\lambda \approx 0{,}1\,\sigma$.

**Algorithme (Chambolle-Pock).** Avec les variables duales $p$ (attache aux données) et $q$ (gradient), et des pas $\tau = \sigma_d = 1/3$ :

$$
\begin{aligned}
p &\leftarrow \frac{p + \sigma_d\,(A\bar{x} - y)}{1 + \sigma_d} \\
q &\leftarrow \Pi_{\lVert q\rVert \le \lambda}\big(q + \sigma_d\,\nabla \bar{x}\big) \\
x^{+} &= \Pi_{[0,1]}\big(x - \tau\,(A^{\top}p - \operatorname{div} q)\big) \\
\bar{x} &\leftarrow 2x^{+} - x
\end{aligned}
$$

Les pas vérifient $\tau\,\sigma_d\,L^2 \le 1$, car $\lVert A\rVert \le 1$ et $\lVert\nabla\rVert^2 \le 8$.

---

## 3. Estimer le noyau quand il est inconnu

### Bougé rectiligne : le cepstre

Le spectre d'un flou rectiligne de longueur $L$ est un sinus cardinal :

$$|K(\omega)| \propto \left|\frac{\sin(\pi L \omega)}{\pi L \omega}\right|$$

Il s'annule régulièrement : le spectre de la photo floue montre des bandes sombres perpendiculaires au mouvement. Le **cepstre**

$$C = \mathcal{F}^{-1}\big(\log |Y|\big)$$

transforme cette périodicité en un pic **négatif** à la distance $\pm L$ de l'origine, dans la direction du mouvement. On cherche le minimum de $C$ dans la couronne $3 \le r \le L_{\max}$. Sa position $(\Delta x, \Delta y)$ donne :

$$L = \sqrt{\Delta x^2 + \Delta y^2}, \qquad \theta = \operatorname{atan2}(-\Delta y,\ \Delta x) \bmod 180^\circ$$

> **En clair :** les bandes sombres sont comme des vagues régulières. Le cepstre mesure leur espacement (la longueur du bougé) et leur orientation (l'angle). Une texture déjà périodique, comme un mur de briques, crée ses propres vagues et trompe la méthode.

### Noyau quelconque : estimation aveugle

Schéma de Xu, Zheng & Jia (2013) et Pan et al. (2014). On part d'un noyau initial en croix, puis on alterne deux étapes, **du grossier au fin**.

#### Étape 1 : image latente avec prior $L_0$

$$\min_x\ \lVert k*x-y\rVert^2 + \lambda\,\lVert\nabla x\rVert_0$$

où $\lVert\nabla x\rVert_0$ compte le nombre de pixels dont le gradient est non nul. On introduit un gradient auxiliaire $g$ (décomposition quadratique) :

$$\min_{x,g}\ \lVert k*x-y\rVert^2 + \beta\,\lVert\nabla x-g\rVert^2 + \lambda\,\lVert g\rVert_0$$

À $\beta$ fixé, les deux sous-problèmes se résolvent exactement :

$$
g = \begin{cases} \nabla x & \text{si } \lVert\nabla x\rVert^2 \ge \lambda/\beta \\ 0 & \text{sinon} \end{cases}
\qquad\qquad
\hat X = \frac{\overline{K}\,Y + \beta\,\mathcal{F}(\nabla^{\top} g)}{|K|^2 + \beta\,|D|^2}
$$

où $D$ est la transformée de Fourier de l'opérateur gradient. On double $\beta$ à chaque tour jusqu'à une valeur très grande, ce qui force $g \to \nabla x$.

> **En clair :** on met à zéro tous les petits gradients (textures, bruit) et on ne garde que les contours francs. Ce sont eux qui renseignent le mieux sur la forme du flou.

#### Étape 2 : le noyau

$$\min_k\ \sum_{d\in\{x,y\}} \lVert k*\partial_d x-\partial_d y\rVert^2 + \gamma\,\lVert k\rVert^2$$

Ce problème a une solution en forme fermée dans Fourier :

$$\hat K = \frac{\sum_d \overline{\mathcal{F}(\partial_d x)}\ \mathcal{F}(\partial_d y)}{\sum_d |\mathcal{F}(\partial_d x)|^2 + \gamma}$$

On travaille sur les gradients plutôt que sur les intensités, car le problème est mieux conditionné. Le résultat est ensuite ramené à un noyau plausible :

1. positivité : $k \leftarrow \max(k, 0)$
2. seuil : mise à zéro des valeurs inférieures à 5 % du maximum
3. suppression des petites composantes connexes isolées (masse inférieure à 10 %)
4. normalisation : $\sum k = 1$
5. recentrage sur le barycentre

#### Du grossier au fin

On construit une pyramide d'images réduites par facteurs successifs de $\sqrt{2}$. Le nombre d'échelles est

$$n_s = \left\lfloor \frac{\ln(k_s/5)}{\ln\sqrt{2}} \right\rfloor + 1$$

où $k_s$ est la taille maximale du noyau. À chaque échelle, on alterne 10 fois les deux étapes, avec $\lambda \leftarrow \max(\lambda/1{,}1,\ 10^{-4})$ pour laisser entrer progressivement plus de détails. Le noyau obtenu est agrandi par interpolation bilinéaire pour initialiser l'échelle suivante.

> **En clair :** sur une petite image, le flou ne fait que quelques pixels et la forme globale du noyau est facile à saisir. On affine ensuite sur des images de plus en plus grandes. Cela évite de se coincer dans une mauvaise solution.

---

## 4. Mesurer la qualité

**PSNR** (images dans $[0,1]$) : plus il est élevé, plus l'image est proche de l'original.

$$\mathrm{PSNR} = 10\log_{10}\frac{1}{\mathrm{MSE}}, \qquad \mathrm{MSE}=\frac1N\sum_i (x_i-\hat x_i)^2$$

**SSIM** : compare luminance, contraste et structure localement, avec une fenêtre gaussienne $11\times 11$ ($\sigma = 1{,}5$).

$$\mathrm{SSIM}=\frac{(2\mu_x\mu_y+C_1)(2\sigma_{xy}+C_2)}{(\mu_x^2+\mu_y^2+C_1)(\sigma_x^2+\sigma_y^2+C_2)}$$

**PSNR recalé.** En aveugle, le noyau n'est connu qu'à une translation près : l'image restaurée peut être décalée de quelques pixels. On retient donc le meilleur PSNR sur les décalages entiers de $\pm 8$ pixels.

**Similarité de noyaux** : corrélation croisée maximale normalisée, avec $s = 1$ pour deux noyaux identiques à une translation près.

$$s(k, k_t) = \frac{\max_\tau\,(k \star k_t)(\tau)}{\lVert k\rVert_2\,\lVert k_t\rVert_2} \in [0,1]$$

---

## 5. Résultats de référence

Bruit gaussien $\sigma = 0{,}01$, images $256 \times 256$.

### Noyau connu (PSNR moyen sur 4 images, en dB)

| Noyau | Flou | Wiener | Richardson-Lucy (40 it.) | TV ($\lambda = 10^{-3}$) |
|---|---|---|---|---|
| K1 linéaire 15 px, 30° | 21,40 | 24,54 | 25,51 | **27,10** |
| K2 linéaire 25 px, 120° | 19,94 | 22,88 | 24,01 | **25,27** |
| K3 tremblement (boucle) | 18,71 | 24,66 | 25,98 | **28,65** |
| K4 tremblement | 18,97 | 26,73 | 27,52 | **30,90** |

La TV est la meilleure partout. Le SSIM moyen passe de 0,49 (flou) à 0,90 (TV) pour K4.

### Cepstre

- K1 (15 px, 30°) : estimé à 15,3 px / 31,6°
- K2 (25 px, 120°) : estimé à 23,7 à 24,2 px / 117,6° à 119,7°
- Sur l'image de briques, la périodicité de la texture trompe la méthode (14,3 px / 24,8° pour K1)

### Estimation aveugle

TV avec le noyau estimé, contre TV avec le vrai noyau (dB), sur 3 images :

| Noyau | Similarité | TV noyau estimé | TV vrai noyau |
|---|---|---|---|
| K1 | 0,92 / 0,91 / 0,93 | 25,1 / 24,0 / 26,2 | 26,0 / 25,5 / 27,6 |
| K2 | 0,90 / 0,90 / 0,90 | 22,7 / 22,1 / 24,2 | 24,2 / 23,9 / 25,1 |
| K3 | 0,85 / 0,79 / 0,71 | 26,0 / 24,6 / 24,8 | 27,5 / 27,3 / 29,0 |
| K4 | 0,91 / 0,91 / 0,88 | 29,1 / 28,1 / 28,2 | 30,9 / 29,2 / 29,8 |

L'estimation aveugle récupère l'essentiel du gain : de 18 dB (flou) à 28-29 dB pour K4, à 1 ou 2 dB de la déconvolution avec le vrai noyau. Elle échoue sur l'image de briques (similarité 0,3 à 0,6) : une texture répétitive sans contour dominant ne contraint pas le noyau.

---

## Limites et pistes

- Le flou est supposé uniforme sur toute l'image : pas de rotation de l'appareil, pas d'objet en mouvement.
- Les PSNR des méthodes aveugles sont mesurés au meilleur recalage.
- Pistes : prior sur les intensités sombres, noyaux variables dans l'image, réseaux de déconvolution déroulés.

