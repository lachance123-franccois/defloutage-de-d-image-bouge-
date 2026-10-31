# Défloutage de bougé

Restaurer une photo devenue floue parce que l'appareil a bougé pendant la prise de vue.

## Le problème

Quand l'appareil bouge pendant la pose, chaque point de la scène laisse une traînée au lieu d'un point net. Chaque pixel de la photo floue est en fait un mélange de ses voisins, le long de la trajectoire de l'appareil.

Cette trajectoire s'appelle le **noyau de flou** : c'est une petite image qui dessine le chemin suivi par l'appareil. À cela s'ajoute un léger bruit, toujours présent sur une photo.

Le projet répond à deux questions :

1. Si on connaît la trajectoire, comment retrouver l'image nette ?
2. Si on ne la connaît pas, comment la deviner à partir de la seule photo floue ?

## Simuler un flou réaliste

**Bougé en ligne droite.** Un simple segment, défini par sa longueur et son angle. Il est dessiné avec soin pour éviter l'effet d'escalier.

**Tremblement aléatoire.** Une courbe lisse et sinueuse, parcourue à vitesse variable. Là où l'appareil ralentit, il reste plus longtemps, donc le noyau y est plus lumineux.

**Les bords de l'image.** Un calcul par transformée de Fourier suppose que l'image se répète comme un carrelage, ce qui fait « déborder » le bord droit sur le bord gauche. Une vraie photo ne fonctionne pas ainsi : la scène continue au-delà du cadre. On prolonge donc l'image en répétant ses pixels de bord, ce qui supprime les défauts parasites sur les contours.

## Retrouver l'image quand le flou est connu

Trois méthodes sont comparées.

**Wiener.** L'idée naturelle est de « défaire » le flou par une division, mais elle amplifie énormément le bruit là où le flou a effacé l'information. Wiener ajoute donc un garde-fou : on ne restaure que ce qui est assez fort pour être fiable. Un réglage dose la prudence.

**Richardson-Lucy.** On part d'une image, on la floute, puis on la compare à la photo réelle. L'écart sert à corriger l'image, et on recommence. Comme les corrections sont multiplicatives, l'image reste toujours positive. Après trop de tours, la méthode se met à reproduire le bruit : la qualité remonte puis redescend. Le nombre d'itérations (une quarantaine ici) est donc un vrai réglage.

**Variation totale.** On cherche l'image qui ressemble à la photo observée tout en étant la plus « simple » possible. Une photo naturelle est faite de zones presque uniformes séparées par des contours nets. La méthode pénalise les petites variations partout, qui sont typiquement du bruit, mais tolère quelques contours francs. Un réglage dose le compromis entre fidélité aux données et lissage. C'est la méthode la plus performante des trois.

## Deviner le flou quand il est inconnu

### Bougé en ligne droite : l'analyse du spectre

Un flou rectiligne laisse une signature dans le spectre de l'image : des bandes sombres régulières, perpendiculaires au mouvement. On peut les voir comme des vagues. Leur espacement donne la longueur du bougé, et leur orientation donne l'angle.

Pour lire cet espacement, on utilise le cepstre, c'est-à-dire le spectre du spectre : il transforme ces vagues régulières en un pic bien net, dont la position indique directement longueur et angle.

Limite : une texture déjà périodique, comme un mur de briques, crée ses propres vagues et trompe la méthode.

### Trajectoire quelconque : estimation aveugle

Pour un tremblement de forme libre, on alterne deux étapes, en répétant le tout du grossier au fin.

**Étape 1 : estimer l'image nette.** On ne garde que les contours francs de l'image et on met à zéro tout le reste (textures fines, bruit). Ce sont les contours qui renseignent le mieux sur la forme du flou : un contour net devenu une traînée montre directement le chemin de l'appareil.

**Étape 2 : estimer la trajectoire.** Connaissant l'image nette approximative, on déduit la trajectoire qui explique le mieux la photo floue. On la nettoie ensuite pour qu'elle ressemble à une vraie trajectoire :
- on supprime les valeurs négatives ;
- on élimine les valeurs trop faibles ;
- on retire les petits fragments isolés ;
- on la normalise ;
- on la recentre.

**Du grossier au fin.** On commence sur une petite version de l'image, où le flou ne fait que quelques pixels et où sa forme globale est facile à saisir. On agrandit ensuite progressivement, en affinant le noyau à chaque taille. Cela évite de s'enfermer dans une mauvaise solution. Au fil des tours, on devient aussi moins sévère sur les contours retenus, pour laisser entrer plus de détails.

## Mesurer la qualité

- **Le PSNR** mesure l'écart entre l'image restaurée et l'originale. Plus il est élevé, mieux c'est.
- **Le SSIM** compare la structure des images (luminosité, contraste, formes), de façon plus proche de la perception humaine.
- **Le PSNR recalé** sert pour les méthodes aveugles : la trajectoire n'étant connue qu'à un décalage près, l'image restaurée peut être légèrement décalée. On cherche donc le meilleur alignement avant de comparer.
- **La similarité de noyaux** compare la trajectoire estimée à la vraie. Elle vaut 1 si elles sont identiques, 0 si elles n'ont aucun rapport.

## Ce qu'on observe

Les tests utilisent 4 images, 4 trajectoires (2 bougés rectilignes, 2 tremblements) et un petit bruit.

- **Flou connu :** la variation totale est la meilleure méthode dans tous les cas, suivie de Richardson-Lucy puis de Wiener. Sur le meilleur cas, la qualité passe d'environ 19 dB (flou) à environ 31 dB.
- **Bougé rectiligne :** l'analyse du spectre retrouve la longueur et l'angle avec une bonne précision, par exemple 15,3 px pour 15 px réels. Elle est trompée par l'image de briques.
- **Flou inconnu :** l'estimation aveugle récupère l'essentiel du gain. Le résultat n'est qu'à 1 ou 2 dB de celui obtenu avec la vraie trajectoire. Elle échoue sur l'image de briques, où la texture répétitive ne laisse apparaître aucun contour dominant.

## Limites et pistes

- Le flou est supposé identique sur toute l'image : pas de rotation de l'appareil, pas d'objet en mouvement.
- Pistes d'amélioration : exploiter les zones sombres de l'image, gérer un flou variable selon les endroits, utiliser des réseaux de neurones construits à partir de ces algorithmes.

