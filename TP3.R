# GOUJON Paul, LOUIS Stephane
# SY09 - UTC - TP3

donn1 <- read.table("Synth1-40.txt", header=F);
donn2 <- read.table("Synth1-100.txt", header=F);
donn3 <- read.table("Synth1-500.txt", header=F);
donn4 <- read.table("Synth1-1000.txt", header=F);
donn5 <- read.table("Synth2-1000.txt", header=F);
donn_list = list(donn1, donn2, donn3, donn4, donn5);
number_experiments = 20;

source("separ1.R");
source("separ2.R");
source("ceuc.R");
source("kppv.R");
individuals_numbers = c(40, 100, 500, 1000, 1000);
# file number suffixe (for the naming of the result files)
fns = c("1", "1", "1", "1", "2");

# SKIP : we wanted the whole TP in the same file. Skip is used to skip some
# parts...
skip = FALSE;

if (!skip) {
# Q 1.2.1
# For each data frame, we want to get the caracteristics of the populations
# Each class is generated following a multi-dimensional normal law
# => we want to estimate mu (the mean) and epsilon (the variance) as well as pi
# (the proportion of each class)
for (i in 1:length(donn_list)) {
    data = as.data.frame(donn_list[i]);
    mu1 = apply(data[which(data$V3 == 1), 1:2], 2, mean);
    mu2 = apply(data[which(data$V3 == 2), 1:2], 2, mean);
    epsilon1 = cov(data[which(data$V3 == 1), 1:2]);
    epsilon2 = cov(data[which(data$V3 == 2), 1:2]);
    pi1 = dim(data[which(data$V3 == 1),])[1] / dim(data)[1];
    pi2 = dim(data[which(data$V3 == 2),])[1] / dim(data)[1];
    df = data.frame("mu1" = mu1, "mu2" = mu2, "epsilon1" = epsilon1, "epsilon2" = epsilon2, "pi1" = pi1, "pi2" = pi2);
    write.table(df, paste("synth", fns[i],"_", individuals_numbers[i], "_analysis.csv", sep=""));
}
} # skip end
# Q 1.2.2

# The aim is to code a script which does 20 data separation, trains a model with
# the training data and test it on the test ensemble in order to calculate the
# error rate

if (!skip) {
ceuc_train_error = matrix(nrow=length(donn_list), ncol=number_experiments);
ceuc_test_error = matrix(nrow=length(donn_list), ncol=number_experiments);

for (j in 1:length(donn_list)) {
    print(paste("CEUC - Dataset: ", j));
    data = as.data.frame(donn_list[j]);
    for (i in 1:number_experiments) {
        print(i);

        # Data separation
        data.sep = separ1(data[,1:2], data[,3]);

        # EUCLIDIAN CLASSIFIER
        mu = ceuc.app(data.sep$Xapp, data.sep$zapp);
        train = ceuc.val(mu, data.sep$Xapp);
        test = ceuc.val(mu, data.sep$Xtst);
        ceuc_train_error[j, i] = dim(as.matrix(which(train != as.matrix(data.sep$zapp))))[1] / dim(data.sep$Xapp)[1];
        ceuc_test_error[j, i] = dim(as.matrix(which(test != as.matrix(data.sep$ztst))))[1] / dim(data.sep$Xtst)[1];
        write.table(ceuc_train_error, "ceuc_train_error.csv");
        write.table(ceuc_test_error, "ceuc_test_error.csv");
        ceuc_train_error_means = apply(ceuc_train_error, 1, mean);
        ceuc_test_error_means = apply(ceuc_test_error, 1, mean);
        write.table(ceuc_train_error_means, "ceuc_train_error_means.csv");
        write.table(ceuc_test_error_means, "ceuc_test_error_means.csv");
    }
}
} # skip end

# CONFIDENCE INTERVALS
# we want to calculate the CIs for each error estimator
# we use the classical alpha = 0.05
# thus we'll need alphas fractil u(1 - (alpha/2)) <=> u(0.975) <=> 1.96

if (!skip) {
u = 1.9600
for (j in 1:length(donn_list)) {
  ptest = ceuc_test_error_means[j];
  ptrain = ceuc_train_error_means[j];
  inftest = ptest - u * sqrt(ptest * (1 - ptest) / individuals_numbers[j]);
  inftrain = ptrain - u * sqrt(ptrain * (1 - ptrain) / individuals_numbers[j]);
  suptest = ptest + u * sqrt(ptest * (1 - ptest) / individuals_numbers[j]);
  suptrain = ptrain + u * sqrt(ptrain * (1 - ptrain) / individuals_numbers[j]);
  ICtest = c(inftest, suptest);
  ICtrain = c(inftrain, suptrain);
  write.table(ICtest, paste("CEUC", fns[j],"_ICtest_", individuals_numbers[j], ".csv", sep=""));
  write.table(ICtrain, paste("CEUC", fns[j],"_ICtrain_", individuals_numbers[j], ".csv", sep=""));
}
} # skip end

# Q 1.2.3
if (!skip) {
nppv = 1:10;
data = as.data.frame(donn_list[4]);
data.sep = separ2(data[,1:2], data[,3]);
kopt = kppv.tune(data.sep$Xapp, data.sep$zapp, data.sep$Xapp, data.sep$zapp, nppv);
write.table(kopt, "kopt1.2.3.csv");
} # end skip

# Q 1.2.4
# Same as 1.2.2. but with the knn algorithm (kppv)
skip = FALSE;
if (!skip) {

kppv_test_error = matrix(nrow=length(donn_list), ncol=number_experiments);
kppv_train_error = matrix(nrow=length(donn_list), ncol=number_experiments);

nppv = 1:10;
for (j in 1:length(donn_list)) {
    print(paste("KPPV - Dataset:", j));
    data = as.data.frame(donn_list[j]);
    for (i in 1:number_experiments) {
        data.sep = separ2(data[,1:2], data[,3]);
        kopt = kppv.tune(data.sep$Xapp, data.sep$zapp, data.sep$Xval, data.sep$zval, nppv);
        ztst = kppv.app(data.sep$Xapp, data.sep$zapp, kopt, data.sep$Xtst);
        zapp = kppv.app(data.sep$Xapp, data.sep$zapp, kopt, data.sep$Xapp);
        kppv_test_error[j, i] = length(ztst[which(as.matrix(ztst) != as.matrix(data.sep$ztst))]) / dim(data.sep$Xtst)[1];
        kppv_train_error[j, i] = length(zapp[which(as.matrix(zapp) != as.matrix(data.sep$zapp))]) / dim(data.sep$Xapp)[1];
    }
}
kppv_test_error_means = apply(kppv_test_error, 1, mean);
kppv_train_error_means = apply(kppv_train_error, 1, mean);
write.table(kppv_test_error, "kppv_test_error.csv");
write.table(kppv_train_error, "kppv_train_error.csv");
write.table(kppv_test_error_means, "kppv_test_error_means.csv");
write.table(kppv_train_error_means, "kppv_train_error_means.csv");
} # skip end

# CONFIDENCE INTERVALS
# we want to calculate the CIs for each error estimator
# we use the classical alpha = 0.05
# thus we'll need alphas fractil u(1 - (alpha/2)) <=> u(0.975) <=> 1.96

if (!skip) {
u = 1.9600
for (j in 1:length(donn_list)) {
  ptest = kppv_test_error_means[j];
  ptrain = kppv_train_error_means[j];
  inftest = ptest - u * sqrt(ptest * (1 - ptest) / individuals_numbers[j]);
  inftrain = ptrain - u * sqrt(ptrain * (1 - ptrain) / individuals_numbers[j]);
  suptest = ptest + u * sqrt(ptest * (1 - ptest) / individuals_numbers[j]);
  suptrain = ptrain + u * sqrt(ptrain * (1 - ptrain) / individuals_numbers[j]);
  ICtest = c(inftest, suptest);
  ICtrain = c(inftrain, suptrain);
  write.table(ICtest, paste("KPPV", fns[j], "_ICtest_", individuals_numbers[j], ".csv", sep=""));
  write.table(ICtrain, paste("KPPV", fns[j], "_ICtrain_", individuals_numbers[j], ".csv", sep=""));
}
} # skip end

