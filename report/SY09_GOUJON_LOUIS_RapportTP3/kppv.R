# GOUJON Paul, LOUIS Stephane
# SY09 - UTC - TP3
# Exercice KPPV


kppv.app = function(Xapp, zapp, K, Xtst) {
    # determines a classification of the individual-variable table passed as
    # argument.
    # IN : 
    # - Xapp : training individuals
    # - zapp : labels of the training individuals
    # - K : numbers of neighbours taken into account
    # - Xtst : test individuals to classificate
    # OUT :
    # - ztst : a matrix n x 1 containing the classification of the n individuals
    # in the Xtst matrix of individiduals

    # n is inidividuals number to test
    n = dim(as.matrix(Xtst))[1];
    # p is individual's caracteristics number
    p = dim(as.matrix(Xapp))[2];
    # nbclass is the number of different classes
    nbclass = length(unique(as.character(zapp)));
    # ztst is the returned vector of classes
    ztst = matrix(nrow=n, ncol=1);
    
    # for individu in Xtst
    for (i in 1:n) {
        # trouver les K plus proches individus de Xapp
        zknn = zapp[order(dist2(Xapp, Xtst[i,]), decreasing=F)][1:K];    
        # trouver la classe la plus représentée
        scores_classes = as.data.frame(table(zknn));
        max_score = max(scores_classes$Freq);
        new_k = K;
        # tant que plusieurs fois score max, refaire en prenant moins de vinzinhos
        while (dim(scores_classes[which(scores_classes$Freq == max_score),])[1] > 1) {
            new_k = new_k + 1;
            zknn = zapp[order(dist2(Xapp, Xtst[i,]), decreasing=F)][1:new_k];    
            scores_classes = as.data.frame(table(zknn));
            max_score = max(scores_classes$Freq);
        }
        # l'attribuer
        ztst[i,] = as.numeric(as.character(scores_classes[which.max(scores_classes$Freq),]$zknn));
    } 
    return (ztst);
}


kppv.tune = function(Xapp, zapp, Xval, zval, nppv) {
    # determines the optimal neighbours number to take into account
    # to determine the individual's class
    # IN:
    # - Xapp : the training individuals matrix
    # - zapp : the labels of the training individuals
    # - Xval : the validation individuals matrix (used to test the optimal value
    # of K)
    # - zval : the validation individual labels
    # - nppv : the different values of K to test
    # OUT:
    # - Kopt : value of the optimal K

    # need library "flexclust" for dist2
    # we use this lib instead of the distXY function because we had some
    # problems using distXY at first because we were not passing arguments
    # casted into the right type (we figured this out afterward)
    library(flexclust);

    # K values matrix
    Kopt = matrix(nrow=dim(as.matrix(nppv))[1], ncol=1);
    # pour chaque valeur de K
    for (i in 1:length(nppv)) {
        # pour chaque individu de X val
        ztst = kppv.app(Xapp, zapp, nppv[i], Xval); 
        # calcul du taux d'erreur
        Kopt[i,] = (length(ztst[which(as.matrix(ztst) != as.matrix(zval)),]) / length(zval) * 100); 
    }
    
    return (nppv[which.min(Kopt)]);
}
