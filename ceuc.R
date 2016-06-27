# Paul GOUJON, Stéphane LOUIS
# UTC - SY09 - TP3
# Exercice 1

ceuc.app = function(Xapp, zapp) {
    # fonction d'apprentissage des paramètres du classifieur euclidien
    # IN :
    # - Xapp : tableau individus - variables Xapp (dimension napp x p)
    # correspondant aux individus d'apprentissage -> matrix
    # - zapp : vecteur de longueur napp des étiquettes associées à chacun des
    # individus -> factor
    # OUT :
    # Retourne les paramètres estimés du classifieur euclidien, sous la forme
    # d'une matrice mu de dimension g x p

    # calcul x1 et x2
    # avec x1 et x2 les centres des nuages des classes 1 et 2
    k = length(unique(as.character(zapp))); 
    mu = matrix(nrow=k, ncol=dim(Xapp)[2]);
    for (i in 1:k) {
        classi = sort(unique(as.character(zapp)))[i]; # class is level(i)
        # individuals are individuals which class IS level(i)
        individualsi = Xapp[which(as.character(zapp) == classi),];
        # x(i) is mean by column of the individuals
        mu[i,] = apply(individualsi, 2, mean);
    }
    return (mu);
}

ceuc.val = function(mu, Xtst) {
    # fonction de classement du tableau individu-variable de test de dimension
    # ntst x p.
    # IN : 
    # - mu : matrice de paramètres estimés
    # - Xtst : ensemble à évaluer
    # OUT : 
    # - vecteur de longueur ntst d'étiquettes prédites.

    # Note : s'aider de la fonction distXY qui calcule les distances
    # euclidiennes au carré entre les individus de deux ensembles X et Y, et de
    # la fonction which.min qui détermine l'indice de l'élément minimal d'un
    # vecteur.
    n = dim(Xtst)[1];
    Ntst = matrix(nrow=n, ncol=1);
    for (i in 1:n) {
        decision = t(mu[2,] - mu[1,]) %*% t(as.matrix(Xtst[i,] - ((mu[1,] + mu[2,])/2)));
        if (decision <= 0) {
            Ntst[i,1] = 1;
        } else {
            Ntst[i,1] = 2;
        }
    }
    return (Ntst);
}

