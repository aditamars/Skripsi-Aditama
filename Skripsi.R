##############################################################################
# PERAMALAN NILAI TUKAR PETANI INDONESIA
# SARIMA, NNAR, DAN HYBRID SARIMA-NNAR
# Estimasi manual menggunakan base R
##############################################################################

# ============================================================================
# 1. DATA DAN PEMBAGIAN TRAINING-TESTING
# ============================================================================

ntp <- c(
  100.69,100.59,98.79,99.02,100.15,100.64,101.71,102.00,101.69,99.20,98.36,98.99,
  98.30,98.77,98.78,99.26,99.41,99.56,99.82,100.24,100.90,100.79,101.13,101.20,
  101.19,101.09,101.20,101.15,101.16,101.39,101.77,101.82,102.19,102.61,102.89,102.75,
  103.01,103.33,103.32,103.91,104.50,104.79,104.87,105.11,105.17,105.51,105.64,105.75,
  105.73,105.10,104.68,104.71,104.77,104.88,104.96,105.26,105.41,105.76,105.72,105.87,
  105.67,105.19,104.53,104.55,104.95,105.28,104.58,104.32,104.56,105.30,105.15,105.30,
  101.95,101.79,101.86,101.80,101.88,101.98,102.12,102.06,102.36,102.87,102.37,101.32,
  101.86,102.19,101.53,100.14,100.02,100.52,100.97,101.28,102.33,102.46,102.95,102.83,
  102.55,102.23,101.32,101.22,101.55,101.47,101.39,101.56,102.02,101.71,101.31,101.49,
  100.91,100.33,99.95,100.01,100.15,100.53,100.65,101.60,102.22,102.78,103.07,103.06,
  102.92,102.33,101.94,101.61,101.99,102.04,101.66,102.56,103.17,103.02,103.12,103.16,
  103.32,102.94,102.73,102.23,102.61,102.33,102.63,103.22,103.88,104.04,104.10,104.46,
  104.16,103.35,102.09,100.32,99.47,99.60,100.09,100.65,101.66,102.25,102.86,103.25,
  103.26,103.10,103.29,102.93,103.39,103.59,103.48,104.68,105.68,106.67,107.18,108.34,
  108.67,108.83,109.29,108.46,105.41,105.96,104.25,106.31,106.82,107.27,107.81,109.00,
  109.84,110.53,110.85,110.58,110.20,110.41,110.64,111.85,114.14,115.78,116.73,117.76,
  118.27,120.97,119.39,116.79,116.71,118.77,119.61,119.85,120.30,120.70,121.29,122.78,
  123.68,123.45,123.72,121.06,121.15,121.72,122.64,123.57,124.36,124.33,124.05,125.35
)

tanggal <- seq(as.Date("2008-01-01"), by="month", length.out=length(ntp))
S <- 12
cat("Jumlah observasi:", length(ntp), "\n")

# ============================================================================
# 2. STATISTIK DESKRIPTIF DAN PEMBAGIAN DATA
# ============================================================================

statistik_deskriptif <- function(y) {
  m <- mean(y); s <- sd(y)
  data.frame(
    n=length(y), mean=m, sd=s, min=min(y), max=max(y),
    median=median(y), skewness=mean((y-m)^3)/s^3,
    kurtosis=mean((y-m)^4)/s^4
  )
}

split_data <- function(y, prop_train=0.90) {
  n <- length(y); n_train <- floor(n*prop_train)
  list(
    train=y[1:n_train],
    test=y[(n_train+1):n],
    n_train=n_train,
    n_test=n-n_train
  )
}

print(statistik_deskriptif(ntp))

plot(tanggal, ntp, type="l", lwd=1.5,
     xlab="Waktu (Bulan)", ylab="Nilai Tukar Petani (NTP)",
     main="NTP Nasional Januari 2008 - Desember 2025")
grid()

sp <- split_data(ntp)
y_train <- sp$train
y_test <- sp$test

cat("Training:", sp$n_train, "| Testing:", sp$n_test, "\n")

# ============================================================================
# 3. UJI POLA MUSIMAN - PERIODOGRAM FISHER'S KAPPA
# ============================================================================

uji_periodogram <- function(y, S=12, alpha=0.05) {
  pgram <- spec.pgram(y, log="no", plot=FALSE)
  I <- pgram$spec
  freq <- pgram$freq
  n <- length(y)
  N <- if(n %% 2 == 0) n/2-1 else (n-1)/2
  
  I_max <- max(I)
  idx_max <- which.max(I)
  T_stat <- I_max/sum(I)
  g_alpha <- 1-(alpha/N)^(1/(N-1))
  
  freq_dominan <- freq[idx_max]
  periode_dominan <- 1/freq_dominan
  
  keputusan <- if(T_stat > g_alpha)
    "Terdapat pola musiman (Tolak H0)"
  else
    "Tidak terdapat pola musiman (Gagal tolak H0)"
  
  freq_musiman <- 1/S
  ymax <- I_max*1.20
  
  plot(freq, I, type="o", pch=16, cex=.6, lwd=1.5,
       col="steelblue4", xlab="Frekuensi (ω)", ylab="Spektrum",
       main="Periodogram Spektral NTP Nasional",
       ylim=c(0,ymax), xaxt="n")
  
  axis(1, at=seq(0,max(freq),by=.05),
       labels=sprintf("%.2f",seq(0,max(freq),by=.05)))
  grid(nx=NA, ny=NULL, col="gray90", lty=1)
  abline(v=freq_musiman, col="red3", lty=2, lwd=2)
  points(freq_dominan, I_max, pch=19, cex=1.2, col="red3")
  abline(h=I_max, col="gray60", lty=3)
  
  text(freq_musiman, ymax*.93,
       paste0("Periode Musiman:\n",S," bulan"),
       col="red3", font=2, cex=.9)
  text(freq_dominan, I_max, paste0("Maks = ",round(I_max,3)),
       pos=4, col="red3", font=2, cex=.8)
  
  mtext(paste0("Frekuensi musiman = 1/",S,
               " = ",round(freq_musiman,4)),
        side=3, line=.2, cex=.8, col="gray30")
  
  list(freq=freq, periodogram=I, I_max=I_max, T=T_stat,
       N=N, g_alpha=g_alpha, lag_dominan=idx_max,
       freq_dominan=freq_dominan, periode_dominan=periode_dominan,
       keputusan=keputusan)
}

hasil_periodogram <- uji_periodogram(y_train, S=12, alpha=.05)

cat("\n=== UJI POLA MUSIMAN ===\n")
cat("Periodogram maksimum :",round(hasil_periodogram$I_max,4),"\n")
cat("Statistik T           :",round(hasil_periodogram$T,4),"\n")
cat("Nilai kritis          :",round(hasil_periodogram$g_alpha,4),"\n")
cat("Frekuensi dominan     :",round(hasil_periodogram$freq_dominan,4),"\n")
cat("Periode dominan       :",round(hasil_periodogram$periode_dominan,2),"bulan\n")
cat("Keputusan             :",hasil_periodogram$keputusan,"\n")

# ============================================================================
# 4. TRANSFORMASI BOX-COX BERTAHAP
# ============================================================================

transformasi_boxcox <- function(y, lambda) {
  if(any(y <= 0)) stop("Data harus bernilai positif.")
  if(abs(lambda)<1e-8) return(log(y))
  (y^lambda-1)/lambda
}

boxcox_loglik <- function(lambda, y) {
  n <- length(y)
  if(any(y <= 0)) return(-Inf)
  
  z <- if(abs(lambda)<1e-8) log(y) else (y^lambda-1)/lambda
  JKS <- sum((z-mean(z))^2)
  
  if(!is.finite(JKS) || JKS <= 0) return(-Inf)
  
  -(n/2)*log(JKS/n)+(lambda-1)*sum(log(y))
}

estimasi_lambda_boxcox <- function(y, interval=c(-2,2),
                                   jumlah_grid=4001) {
  if(any(y <= 0)) stop("Data harus bernilai positif.")
  
  lambda_grid <- seq(interval[1],interval[2],length.out=jumlah_grid)
  loglik <- sapply(lambda_grid,boxcox_loglik,y=y)
  i <- which.max(loglik)
  
  list(lambda=lambda_grid[i],loglik=loglik[i],
       lambda_grid=lambda_grid,loglik_grid=loglik)
}

buat_positif <- function(y) {
  minimum <- min(y)
  shift <- if(minimum <= 0) abs(minimum)+1e-6 else 0
  
  list(data=y+shift,shift=shift)
}

boxcox_bertahap <- function(y, maksimum_transformasi=5,
                            interval_lambda=c(-2,2),
                            toleransi_lambda=.10) {
  
  hasil <- data.frame(Transformasi=integer(),Lambda=numeric(),
                      Mean=numeric(),SD=numeric(),CV=numeric(),
                      Shift=numeric())
  data_transformasi <- list()
  y_sekarang <- y
  
  h <- estimasi_lambda_boxcox(y_sekarang,interval_lambda)
  lambda <- h$lambda
  
  hasil <- rbind(hasil,data.frame(
    Transformasi=0,Lambda=lambda,Mean=mean(y_sekarang),
    SD=sd(y_sekarang),CV=sd(y_sekarang)/abs(mean(y_sekarang)),
    Shift=0))
  data_transformasi[[1]] <- y_sekarang
  
  if(abs(lambda-1)<=toleransi_lambda)
    return(list(ringkasan=hasil,data=data_transformasi,
                data_akhir=y_sekarang))
  
  for(i in 1:maksimum_transformasi) {
    positif <- buat_positif(y_sekarang)
    y_input <- positif$data
    shift <- positif$shift
    
    h <- estimasi_lambda_boxcox(y_input,interval_lambda)
    lambda <- h$lambda
    y_baru <- transformasi_boxcox(y_input,lambda)
    
    mean_baru <- mean(y_baru)
    sd_baru <- sd(y_baru)
    cv_baru <- sd_baru/abs(mean_baru)
    
    hasil <- rbind(hasil,data.frame(
      Transformasi=i,Lambda=lambda,Mean=mean_baru,
      SD=sd_baru,CV=cv_baru,Shift=shift))
    data_transformasi[[i+1]] <- y_baru
    y_sekarang <- y_baru
    
    cat("\nTransformasi ke-",i,
        "\nLambda:",round(lambda,6),
        "\nMean:",round(mean_baru,6),
        "\nSD:",round(sd_baru,6),
        "\nCV:",round(cv_baru,6))
    
    if(shift > 0) cat("\nShift:",round(shift,6))
    
    if(abs(lambda-1)<=toleransi_lambda) {
      cat("\nStatus: Lambda sudah mendekati 1.\n")
      break
    } else cat("\nStatus: Lambda belum mendekati 1.\n")
  }
  
  list(ringkasan=hasil,data=data_transformasi,data_akhir=y_sekarang)
}

hasil_boxcox <- boxcox_bertahap(
  y_train,maksimum_transformasi=5,
  interval_lambda=c(-2,2),toleransi_lambda=.10
)

cat("\n\n=== HASIL TRANSFORMASI BOX-COX ===\n")
print(hasil_boxcox$ringkasan,row.names=FALSE)

# ============================================================================
# 5. STASIONERITAS DALAM RATAAN
# ============================================================================

y_bc <- hasil_boxcox$data_akhir
S <- 12

acf_manual <- function(y,max_lag=48) {
  n <- length(y); y_bar <- mean(y)
  gamma0 <- sum((y-y_bar)^2)/n
  rho <- numeric(max_lag)
  
  for(k in 1:max_lag) {
    gamma_k <- sum((y[(k+1):n]-y_bar)*
                     (y[1:(n-k)]-y_bar))/n
    rho[k] <- gamma_k/gamma0
  }
  rho
}

pacf_manual <- function(y,max_lag=48) {
  rho <- c(1,acf_manual(y,max_lag))
  phi <- matrix(0,max_lag,max_lag)
  pacf <- numeric(max_lag)
  
  phi[1,1] <- rho[2]
  pacf[1] <- rho[2]
  
  for(k in 2:max_lag) {
    phi[k,k] <- (rho[k+1]-
                   sum(phi[k-1,1:(k-1)]*rho[k:2]))/
      (1-sum(phi[k-1,1:(k-1)]*rho[2:k]))
    
    for(j in 1:(k-1))
      phi[k,j] <- phi[k-1,j]-phi[k,k]*phi[k-1,k-j]
    
    pacf[k] <- phi[k,k]
  }
  pacf
}

batas_signifikansi <- function(n) 1.96/sqrt(n)

# ACF dan PACF sebelum differencing
acf_awal <- acf_manual(y_bc)
pacf_awal <- pacf_manual(y_bc)
batas_awal <- batas_signifikansi(length(y_bc))

par(mfrow=c(1,2))
plot(acf_awal,type="h",lwd=2,main="(a) ACF Sebelum Differencing",
     xlab="Lag",ylab="ACF",ylim=c(-1,1))
abline(h=c(-batas_awal,batas_awal),lty=2)

plot(pacf_awal,type="h",lwd=2,main="(b) PACF Sebelum Differencing",
     xlab="Lag",ylab="PACF",ylim=c(-1,1))
abline(h=c(-batas_awal,batas_awal),lty=2)
par(mfrow=c(1,1))

# Differencing SARIMA
diff_manual <- function(y,lag=1,differences=1) {
  for(i in 1:differences) {
    n <- length(y)
    y <- y[(lag+1):n]-y[1:(n-lag)]
  }
  y
}

d_order <- 1
D_order <- 1

w <- diff_manual(
  diff_manual(y_bc,lag=S,differences=D_order),
  lag=1,differences=d_order
)

# ACF dan PACF setelah differencing
acf_w <- acf_manual(w)
pacf_w <- pacf_manual(w)
batas_w <- batas_signifikansi(length(w))

par(mfrow=c(1,2))
plot(acf_w,type="h",lwd=2,main="(a) ACF Setelah Differencing",
     xlab="Lag",ylab="ACF",ylim=c(-1,1))
abline(h=c(-batas_w,batas_w),lty=2)

plot(pacf_w,type="h",lwd=2,main="(b) PACF Setelah Differencing",
     xlab="Lag",ylab="PACF",ylim=c(-1,1))
abline(h=c(-batas_w,batas_w),lty=2)
par(mfrow=c(1,1))

# ============================================================================
# 6. UJI ADF
# ============================================================================

tabel_DF <- data.frame(
  n=c(25,50,100,250,500,10000),
  crit=c(-2.986,-2.921,-2.891,-2.873,-2.867,-2.863)
)

nilai_kritis_DF <- function(n)
  approx(tabel_DF$n,tabel_DF$crit,xout=n,rule=2)$y

adf_test_manual <- function(y) {
  n <- length(y)
  dY <- diff_manual(y)
  X <- cbind(1,y[-n])
  
  beta <- solve(t(X)%*%X)%*%t(X)%*%dY
  resid <- dY-X%*%beta
  
  sigma2 <- sum(resid^2)/(length(dY)-2)
  var_beta <- sigma2*solve(t(X)%*%X)
  
  psi_hat <- beta[2]
  se_psi <- sqrt(var_beta[2,2])
  tau <- psi_hat/se_psi
  crit <- nilai_kritis_DF(n)
  
  keputusan <- if(tau < crit)
    "Stasioner dalam rataan (Tolak H0)"
  else
    "Belum stasioner dalam rataan (Gagal tolak H0)"
  
  list(psi_hat=psi_hat,se_psi=se_psi,tau=tau,
       nilai_kritis=crit,keputusan=keputusan)
}

adf_awal <- adf_test_manual(y_bc)
adf_akhir <- adf_test_manual(w)

cat("\n=== ADF SEBELUM DIFFERENCING ===\n")
cat("Psi:",round(adf_awal$psi_hat,6),
    "| SE:",round(adf_awal$se_psi,6),
    "| ADF:",round(adf_awal$tau,6),
    "| Kritis:",round(adf_awal$nilai_kritis,6),"\n")
cat("Keputusan:",adf_awal$keputusan,"\n")

cat("\n=== ADF SETELAH DIFFERENCING ===\n")
cat("Psi:",round(adf_akhir$psi_hat,6),
    "| SE:",round(adf_akhir$se_psi,6),
    "| ADF:",round(adf_akhir$tau,6),
    "| Kritis:",round(adf_akhir$nilai_kritis,6),"\n")
cat("Keputusan:",adf_akhir$keputusan,"\n")

# ============================================================================
# 7. ESTIMASI SARIMA SECARA MANUAL
# ============================================================================

poly_mult <- function(a,b) {
  z <- numeric(length(a)+length(b)-1)
  for(i in seq_along(a))
    z[i:(i+length(b)-1)] <- z[i:(i+length(b)-1)]+a[i]*b
  z
}

gabung_koef <- function(nonmusiman,musiman,S) {
  ar <- c(1,-nonmusiman)
  
  seas <- if(length(musiman)==0) 1 else {
    z <- numeric(length(musiman)*S+1)
    z[1] <- 1
    z[seq_along(musiman)*S+1] <- -musiman
    z
  }
  
  -poly_mult(ar,seas)[-1]
}

sarima_residual <- function(par,w,p,q,P,Q,S) {
  k <- c(p,q,P,Q)
  idx <- c(0,cumsum(k))
  
  phi <- if(p>0) par[(idx[1]+1):idx[2]] else numeric(0)
  theta <- if(q>0) par[(idx[2]+1):idx[3]] else numeric(0)
  Phi <- if(P>0) par[(idx[3]+1):idx[4]] else numeric(0)
  Theta <- if(Q>0) par[(idx[4]+1):idx[5]] else numeric(0)
  
  ar <- gabung_koef(phi,Phi,S)
  ma <- gabung_koef(theta,Theta,S)
  e <- numeric(length(w))
  
  for(t in seq_along(w)) {
    na <- min(length(ar),t-1)
    nm <- min(length(ma),t-1)
    
    ar_term <- if(na>0) sum(ar[1:na]*w[t-(1:na)]) else 0
    ma_term <- if(nm>0) sum(ma[1:nm]*e[t-(1:nm)]) else 0
    
    e[t] <- w[t]-ar_term-ma_term
  }
  
  awal <- max(length(ar),length(ma))+1
  e[awal:length(e)]
}

negloglik <- function(par,w,p,q,P,Q,S) {
  e <- sarima_residual(par,w,p,q,P,Q,S)
  sigma2 <- mean(e^2)
  
  if(!is.finite(sigma2) || sigma2<=0) return(1e10)
  
  .5*length(e)*(log(2*pi)+log(sigma2)+1)
}

fit_sarima_manual <- function(y,p,d,q,P,D,Q,S) {
  
  # Differencing
  w <- diff_manual(diff_manual(y,lag=S,differences=D),
                   lag=1,differences=d)
  k <- p+q+P+Q
  
  fungsi_objektif <- function(par)
    negloglik(par,w,p,q,P,Q,S)
  
  # Estimasi parameter
  opt <- optim(
    par=rep(.1,k),fn=fungsi_objektif,
    method="L-BFGS-B",
    lower=rep(-.98,k),upper=rep(.98,k),
    hessian=TRUE
  )
  
  e <- sarima_residual(opt$par,w,p,q,P,Q,S)
  
  V <- tryCatch(solve(opt$hessian),
                error=function(err) matrix(NA,k,k))
  se <- sqrt(pmax(diag(V),0))
  z <- opt$par/se
  p_value <- 2*pnorm(-abs(z))
  
  loglik <- -opt$value
  n <- length(e)
  AIC <- -2*loglik+2*(k+1)
  BIC <- -2*loglik+log(n)*(k+1)
  
  list(param=opt$par,se=se,z=z,p_value=p_value,
       residual=e,w=w,loglik=loglik,AIC=AIC,BIC=BIC,
       p=p,d=d,q=q,P=P,D=D,Q=Q,S=S,y_asli=y)
}

# ============================================================================
# 8. KANDIDAT MODEL SARIMA
# ============================================================================

kandidat_model <- list(
  c(p=0, q=0, P=1, Q=0),
  c(p=0, q=0, P=0, Q=1),
  c(p=0, q=0, P=1, Q=1),
  c(p=0, q=1, P=0, Q=1),
  c(p=1, q=0, P=1, Q=0),
  c(p=1, q=1, P=0, Q=1),
  c(p=1, q=1, P=1, Q=1),
  c(p=2, q=0, P=0, Q=1),
  c(p=2, q=0, P=1, Q=0),
  c(p=2, q=1, P=0, Q=1),
  c(p=2, q=1, P=1, Q=0),
  c(p=3, q=0, P=0, Q=1),
  c(p=3, q=0, P=1, Q=0),
  c(p=3, q=1, P=0, Q=1),
  c(p=3, q=1, P=1, Q=0),
  c(p=3, q=1, P=1, Q=1)
)

hasil_kandidat <- lapply(kandidat_model,function(o)
  fit_sarima_manual(
    y_bc,p=o["p"],d=d_order,q=o["q"],
    P=o["P"],D=D_order,Q=o["Q"],S=S
  ))

# ============================================================================
# 9. SIGNIFIKANSI PARAMETER
# ============================================================================

tabel_signifikansi <- do.call(rbind,lapply(seq_along(hasil_kandidat),function(i) {
  f <- hasil_kandidat[[i]]
  o <- kandidat_model[[i]]
  
  data.frame(
    Model=sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]",
                  o["p"],d_order,o["q"],o["P"],D_order,o["Q"],S),
    Parameter=seq_along(f$param),
    Estimasi=f$param,SE=f$se,Z=f$z,
    p_value=f$p_value,
    Signifikan=ifelse(f$p_value<.05,"Ya","Tidak")
  )
}))

print(tabel_signifikansi,row.names=FALSE)

# ============================================================================
# 10. DIAGNOSTIK RESIDUAL
# ============================================================================

ljung_box_manual <- function(e, lag=24, fitdf=0, alpha=.05) {
  e <- e[is.finite(e)]
  n <- length(e)
  rho <- acf(e, lag.max=lag, plot=FALSE)$acf[-1]
  k <- 1:lag
  
  Q <- n*(n+2)*sum(rho^2/(n-k))
  df <- lag-fitdf
  Q_kritis <- qchisq(1-alpha, df)
  p_value <- 1-pchisq(Q, df)
  
  list(
    Q=Q,
    Q_kritis=Q_kritis,
    p_value=p_value,
    keputusan=ifelse(
      p_value>alpha,
      "Residual white noise",
      "Residual tidak white noise"
    )
  )
}

ks_test_manual <- function(e, alpha=.05) {
  e <- e[is.finite(e)]
  n <- length(e)
  
  z <- sort((e-mean(e))/sd(e))
  Fn1 <- (1:n)/n
  Fn2 <- (0:(n-1))/n
  Fz <- pnorm(z)
  
  D <- max(max(Fn1-Fz), max(Fz-Fn2))
  D_alpha <- 1.36/sqrt(n)
  
  lambda <- (sqrt(n)+.12+.11/sqrt(n))*D
  
  p_value <- sum(
    (-1)^(1:100-1)*exp(-2*(1:100)^2*lambda^2)
  )
  p_value <- max(min(2*p_value,1),0)
  
  list(
    D=D,
    D_alpha=D_alpha,
    p_value=p_value,
    keputusan=ifelse(
      p_value>alpha,
      "Residual normal",
      "Residual tidak normal"
    )
  )
}

diagnostik_model <- function(e, fitdf=0, lag=24) {
  lb <- ljung_box_manual(e, lag, fitdf)
  ks <- ks_test_manual(e)
  
  c(
    LB_Q=lb$Q,
    LB_Qkritis=lb$Q_kritis,
    LB_p=lb$p_value,
    KS_D=ks$D,
    KS_Dkritis=ks$D_alpha,
    KS_p=ks$p_value
  )
}

# ---------------------------------------------------------------------------
# Memilih model yang seluruh parameternya signifikan
# ---------------------------------------------------------------------------

model_lolos <- sapply(seq_along(hasil_kandidat), function(i) {
  all(hasil_kandidat[[i]]$p_value < .05)
})

indeks_terpilih <- which(model_lolos)

hasil_signifikan <- hasil_kandidat[indeks_terpilih]
model_signifikan <- kandidat_model[indeks_terpilih]

tabel_diagnostik <- do.call(rbind,
                            lapply(seq_along(hasil_signifikan), function(i) {
                              
                              f <- hasil_signifikan[[i]]
                              o <- model_signifikan[[i]]
                              
                              d <- diagnostik_model(
                                f$residual,
                                length(f$param),
                                24
                              )
                              
                              data.frame(
                                Model=sprintf(
                                  "SARIMA(%d,%d,%d)(%d,%d,%d)[%d]",
                                  o["p"],d_order,o["q"],
                                  o["P"],D_order,o["Q"],S
                                ),
                                Q=round(d["LB_Q"],4),
                                Q_kritis=round(d["LB_Qkritis"],4),
                                LB_p=round(d["LB_p"],6),
                                D_KS=round(d["KS_D"],6),
                                D_alpha_n=round(d["KS_Dkritis"],6),
                                KS_p=round(d["KS_p"],6)
                              )
                            })
)

print(tabel_diagnostik,row.names=FALSE)

# ============================================================================
# 11. AIC DAN BIC
# ============================================================================

tabel_sarima <- do.call(rbind,
                        lapply(seq_along(hasil_signifikan), function(i) {
                          
                          f <- hasil_signifikan[[i]]
                          o <- model_signifikan[[i]]
                          
                          data.frame(
                            Model=sprintf(
                              "SARIMA(%d,%d,%d)(%d,%d,%d)[%d]",
                              o["p"],d_order,o["q"],
                              o["P"],D_order,o["Q"],S
                            ),
                            AIC=round(f$AIC,4),
                            BIC=round(f$BIC,4)
                          )
                        })
)

print(tabel_sarima,row.names=FALSE)

# ============================================================================
# 12. MODEL SARIMA TERBAIK
# ============================================================================

idx <- which.min(tabel_sarima$AIC)

sarima_terbaik <- hasil_signifikan[[idx]]
model_terbaik <- model_signifikan[[idx]]

cat(
  "\nModel terbaik:",
  sprintf(
    "SARIMA(%d,%d,%d)(%d,%d,%d)[%d]",
    model_terbaik["p"],d_order,model_terbaik["q"],
    model_terbaik["P"],D_order,model_terbaik["Q"],S
  ),
  "\nAIC:",round(sarima_terbaik$AIC,4),
  "\nBIC:",round(sarima_terbaik$BIC,4),"\n"
)

# ============================================================================
# 13. PERAMALAN DAN EVALUASI SARIMA TERBAIK
# ============================================================================

# Model terbaik: SARIMA(0,1,0)(0,1,1)[12]

theta <- sarima_terbaik$param[1]
z <- y_bc
h <- length(y_test)

# Double differencing
w <- diff_manual(
  diff_manual(z,lag=12),
  lag=1
)

# Residual
e <- numeric(length(w))

for(t in seq_along(w)) {
  e[t] <- w[t]-
    ifelse(t>12,theta*e[t-12],0)
}

# Peramalan pada skala Box-Cox
z_ext <- c(z,rep(NA,h))

for(i in 1:h) {
  
  t <- length(z)+i
  
  z_ext[t] <-
    z_ext[t-1]+
    z_ext[t-12]-
    z_ext[t-13]+
    ifelse(
      i<=12,
      theta*e[length(e)-12+i],
      0
    )
}

ramalan_bc <- z_ext[
  (length(z)+1):(length(z)+h)
]

# ============================================================================
# 14. INVERS BOX-COX
# ============================================================================

inverse_boxcox <- function(z,lambda) {
  
  if(abs(lambda)<1e-8)
    exp(z)
  else
    (lambda*z+1)^(1/lambda)
}

ramalan_sarima <- ramalan_bc

# Membalik transformasi dari terakhir ke pertama
for(i in (nrow(hasil_boxcox$ringkasan)-1):1) {
  
  lambda <- hasil_boxcox$ringkasan$Lambda[i+1]
  shift <- hasil_boxcox$ringkasan$Shift[i+1]
  
  ramalan_sarima <-
    inverse_boxcox(ramalan_sarima,lambda)-shift
}

# ============================================================================
# 15. EVALUASI MODEL SARIMA TERBAIK
# ============================================================================

error <- y_test-ramalan_sarima

RMSE <- sqrt(mean(error^2))
MAPE <- mean(abs(error/y_test))*100

tanggal_test <-
  tanggal[(sp$n_train+1):length(ntp)]

tabel_peramalan <- data.frame(
  Periode=format(tanggal_test,"%Y-%m"),
  Peramalan=round(ramalan_sarima,3),
  Aktual=round(y_test,3)
)

cat("\n=== HASIL PERAMALAN MODEL SARIMA TERBAIK ===\n")

cat(
  "Model : SARIMA(0,1,0)(0,1,1)[12]\n\n"
)

print(tabel_peramalan,row.names=FALSE)

cat("\nRMSE :",round(RMSE,4),"\n")
cat("MAPE :",round(MAPE,4),"%\n")

# ============================================================================
# NNAR (NEURAL NETWORK AUTOREGRESSION)
# ============================================================================

# ============================================================================
# 1. NORMALISASI DATA
# ============================================================================

# ----------------------------------------------------------------------------
# 1.1 Fungsi Normalisasi Min-Max
# ----------------------------------------------------------------------------

normalisasi_minmax <- function(y, batas = c(0, 1)) {
  y_min <- min(y)
  y_max <- max(y)
  a <- batas[1]
  b <- batas[2]
  
  y_norm <- a + (y - y_min) * (b - a) / (y_max - y_min)
  
  list(
    data = y_norm,
    min = y_min,
    max = y_max,
    batas = batas
  )
}

norm_train <- normalisasi_minmax(y_train)

y_train_norm <- norm_train$data
y_min_train  <- norm_train$min
y_max_train  <- norm_train$max

# Normalisasi data testing menggunakan parameter training
y_test_norm <- (y_test - y_min_train) /
  (y_max_train - y_min_train)


# ----------------------------------------------------------------------------
# 1.2 Lampiran Data Normalisasi
# ----------------------------------------------------------------------------

buat_lampiran_normalisasi <- function(y_norm, jumlah_kolom = 6) {
  
  n <- length(y_norm)
  n_baris <- ceiling(n / jumlah_kolom)
  y_pad <- c(y_norm, rep(NA, n_baris * jumlah_kolom - n))
  
  tabel <- matrix(
    round(y_pad, 4),
    nrow = n_baris,
    ncol = jumlah_kolom
  )
  
  colnames(tabel) <- paste0("Norm_", 1:jumlah_kolom)
  as.data.frame(tabel)
}

lampiran3 <- buat_lampiran_normalisasi(y_train_norm)

cat("\n=== LAMPIRAN 3. DATA NORMALISASI NTP TRAINING ===\n")
print(lampiran3, row.names = FALSE, na.print = "")

write.csv(
  lampiran3,
  "Lampiran3_Normalisasi_Data_Training_NTP.csv",
  row.names = FALSE,
  na = ""
)


# ============================================================================
# 2. PEMBENTUKAN DATA LAG NNAR
# ============================================================================

buat_lag_nnar <- function(y_norm, p, P, S) {
  
  n <- length(y_norm)
  
  lag_nonmusiman <- 1:p
  lag_musiman <- if (P > 0) (1:P) * S else numeric(0)
  lag_semua <- c(lag_nonmusiman, lag_musiman)
  lag_maks <- max(lag_semua)
  
  hasil <- matrix(
    NA,
    nrow = n - lag_maks,
    ncol = length(lag_semua) + 1
  )
  
  for (i in (lag_maks + 1):n) {
    hasil[i - lag_maks, ] <- c(
      y_norm[i],
      y_norm[i - lag_semua]
    )
  }
  
  nama_kolom <- c(
    "Yt",
    paste0("Yt_", lag_nonmusiman),
    if (P > 0) paste0("Yt_", lag_musiman, "(S)") else NULL
  )
  
  df <- data.frame(
    t = (lag_maks + 1):n,
    hasil
  )
  
  colnames(df) <- c("t", nama_kolom)
  
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) round(x, 4) else x
  })
  
  df
}

# Struktur NNAR berdasarkan hasil identifikasi
p_input <- 3
P_input <- 1
S <- 12

lampiran4 <- buat_lag_nnar(
  y_train_norm,
  p = p_input,
  P = P_input,
  S = S
)

cat(
  "\n=== LAMPIRAN 4. DATA PRAPROSES NNAR(",
  p_input, ",", P_input, ")_",
  S, " ===\n",
  sep = ""
)

print(lampiran4, row.names = FALSE)

write.csv(
  lampiran4,
  "Lampiran4_Data_Praproses_NNAR_NTP.csv",
  row.names = FALSE
)


# ============================================================================
# 3. PERANCANGAN JARINGAN NNAR
# ============================================================================

# ----------------------------------------------------------------------------
# 3.1 Fungsi Aktivasi Sigmoid
# ----------------------------------------------------------------------------

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

turunan_sigmoid <- function(y) {
  y * (1 - y)
}


# ----------------------------------------------------------------------------
# 3.2 Inisialisasi Bobot dan Bias
# ----------------------------------------------------------------------------

inisialisasi_bobot <- function(n_input, n_hidden,
                               seed = 123,
                               batas = c(-0.5, 0.5)) {
  
  set.seed(seed)
  
  list(
    v = matrix(
      runif(n_input * n_hidden, batas[1], batas[2]),
      nrow = n_input
    ),
    v0 = runif(n_hidden, batas[1], batas[2]),
    w = matrix(
      runif(n_hidden, batas[1], batas[2]),
      nrow = n_hidden
    ),
    w0 = runif(1, batas[1], batas[2])
  )
}


# ----------------------------------------------------------------------------
# 3.3 Fungsi Pelatihan Backpropagation
# ----------------------------------------------------------------------------

latih_nnar <- function(X, Y, n_hidden, learning_rate,
                       threshold = 0.1,
                       max_epoch = 10000,
                       seed = 123) {
  
  n_data <- nrow(X)
  n_input <- ncol(X)
  
  b <- inisialisasi_bobot(
    n_input,
    n_hidden,
    seed
  )
  
  v <- b$v
  v0 <- b$v0
  w <- b$w
  w0 <- b$w0
  
  epoch <- 0
  sse <- Inf
  
  while (sse > threshold && epoch < max_epoch) {
    
    sse <- 0
    
    for (i in seq_len(n_data)) {
      
      x <- as.numeric(X[i, ])
      t <- Y[i]
      
      # Feed forward
      z_in <- v0 + x %*% v
      z <- sigmoid(z_in)
      
      y_in <- w0 + z %*% w
      y <- sigmoid(y_in)
      
      # Backpropagation
      error <- t - y
      
      delta_output <- error *
        turunan_sigmoid(as.numeric(y))
      
      delta_hidden <- as.numeric(delta_output) *
        as.numeric(w) *
        turunan_sigmoid(as.numeric(z))
      
      # Update bobot
      w <- w +
        learning_rate *
        as.numeric(delta_output) *
        as.numeric(z)
      
      w0 <- w0 +
        learning_rate *
        as.numeric(delta_output)
      
      v <- v +
        learning_rate *
        outer(x, delta_hidden)
      
      v0 <- v0 +
        learning_rate *
        delta_hidden
      
      sse <- sse + as.numeric(error)^2
    }
    
    epoch <- epoch + 1
  }
  
  list(
    v = v,
    v0 = v0,
    w = w,
    w0 = w0,
    epoch = epoch,
    sse_akhir = sse,
    konvergen = sse <= threshold
  )
}


# ----------------------------------------------------------------------------
# 3.4 Fungsi Prediksi
# ----------------------------------------------------------------------------

prediksi_nnar <- function(X, model) {
  
  hasil <- numeric(nrow(X))
  
  for (i in seq_len(nrow(X))) {
    
    x <- as.numeric(X[i, ])
    
    z <- sigmoid(
      model$v0 + x %*% model$v
    )
    
    hasil[i] <- sigmoid(
      model$w0 + z %*% model$w
    )
  }
  
  hasil
}


# ----------------------------------------------------------------------------
# 3.5 Fungsi Evaluasi dan Denormalisasi
# ----------------------------------------------------------------------------

hitung_rmse <- function(aktual, prediksi) {
  sqrt(mean((aktual - prediksi)^2))
}

hitung_mape <- function(aktual, prediksi) {
  mean(abs((aktual - prediksi) / aktual)) * 100
}

denormalisasi_minmax <- function(x, min_data, max_data) {
  x * (max_data - min_data) + min_data
}


# ============================================================================
# 4. GRID SEARCH NNAR
# ============================================================================

# ----------------------------------------------------------------------------
# 4.1 Parameter Grid Search
# ----------------------------------------------------------------------------

neuron_hidden_list <- 2:3

learning_rate_list <- seq(
  0.01, 0.10, by = 0.01
)

threshold_nnar <- 0.1
max_epoch_nnar <- 10000


# ----------------------------------------------------------------------------
# 4.2 Data Input dan Target
# ----------------------------------------------------------------------------

kolom_input <- setdiff(
  colnames(lampiran4),
  c("t", "Yt")
)

X_train <- as.matrix(
  lampiran4[, kolom_input]
)

Y_train <- lampiran4$Yt


# ----------------------------------------------------------------------------
# 4.3 Proses Grid Search
# ----------------------------------------------------------------------------

hasil_grid <- data.frame(
  Neuron_Hidden = integer(),
  Learning_Rate = numeric(),
  RMSE = numeric(),
  MAPE = numeric()
)

for (nh in neuron_hidden_list) {
  
  for (lr in learning_rate_list) {
    
    cat(
      "Neuron Hidden =", nh,
      "| Learning Rate =", lr, "\n"
    )
    
    model <- latih_nnar(
      X = X_train,
      Y = Y_train,
      n_hidden = nh,
      learning_rate = lr,
      threshold = threshold_nnar,
      max_epoch = max_epoch_nnar
    )
    
    pred_norm <- prediksi_nnar(
      X_train,
      model
    )
    
    pred_asli <- denormalisasi_minmax(
      pred_norm,
      y_min_train,
      y_max_train
    )
    
    akt_asli <- denormalisasi_minmax(
      Y_train,
      y_min_train,
      y_max_train
    )
    
    hasil_grid <- rbind(
      hasil_grid,
      data.frame(
        Neuron_Hidden = nh,
        Learning_Rate = lr,
        RMSE = round(
          hitung_rmse(akt_asli, pred_asli),
          4
        ),
        MAPE = round(
          hitung_mape(akt_asli, pred_asli),
          4
        )
      )
    )
  }
}


# ============================================================================
# 5. HASIL TRIAL AND ERROR LEARNING RATE
# ============================================================================

buat_lampiran_lr <- function(df, neuron_list) {
  
  tabel <- data.frame(
    Learning_Rate = learning_rate_list
  )
  
  for (nh in neuron_list) {
    
    sub <- df[df$Neuron_Hidden == nh, ]
    
    tabel[[paste0("RMSE_", nh)]] <- sub$RMSE
    tabel[[paste0("MAPE_", nh)]] <- sub$MAPE
  }
  
  tabel
}

lampiran9 <- buat_lampiran_lr(
  hasil_grid,
  neuron_hidden_list
)

cat("\n=== LAMPIRAN 9. HASIL LEARNING RATE NNAR ===\n")
print(lampiran9, row.names = FALSE)

write.csv(
  lampiran9,
  "Lampiran9_Learning_Rate_NNAR.csv",
  row.names = FALSE
)


# ============================================================================
# 6. PENENTUAN MODEL TERBAIK
# ============================================================================

# Model terbaik berdasarkan RMSE
tabel4_rmse <- do.call(
  rbind,
  lapply(neuron_hidden_list, function(nh) {
    
    sub <- hasil_grid[
      hasil_grid$Neuron_Hidden == nh,
    ]
    
    sub[which.min(sub$RMSE), ]
  })
)

# Model terbaik berdasarkan MAPE
tabel4_mape <- do.call(
  rbind,
  lapply(neuron_hidden_list, function(nh) {
    
    sub <- hasil_grid[
      hasil_grid$Neuron_Hidden == nh,
    ]
    
    sub[which.min(sub$MAPE), ]
  })
)

cat("\n=== TABEL 4A. MODEL TERBAIK BERDASARKAN RMSE ===\n")
print(tabel4_rmse, row.names = FALSE)

cat("\n=== TABEL 4B. MODEL TERBAIK BERDASARKAN MAPE ===\n")
print(tabel4_mape, row.names = FALSE)

write.csv(
  tabel4_rmse,
  "Tabel4_NNAR_Terbaik_RMSE.csv",
  row.names = FALSE
)

write.csv(
  tabel4_mape,
  "Tabel4_NNAR_Terbaik_MAPE.csv",
  row.names = FALSE
)


# ============================================================================
# 7. INISIALISASI BOBOT AWAL MODEL TERBAIK
# ============================================================================

n_input <- ncol(X_train)

idx_terbaik <- which.min(hasil_grid$MAPE)

n_hidden_terbaik <- hasil_grid$Neuron_Hidden[idx_terbaik]
lr_terbaik <- hasil_grid$Learning_Rate[idx_terbaik]

cat("\n=== STRUKTUR JARINGAN TERBAIK ===\n")
cat("Neuron input  :", n_input, "\n")
cat("Neuron hidden :", n_hidden_terbaik, "\n")
cat("Learning rate :", lr_terbaik, "\n")

bobot_awal <- inisialisasi_bobot(
  n_input = n_input,
  n_hidden = n_hidden_terbaik,
  seed = 123
)


# ----------------------------------------------------------------------------
# 7.1 Bobot Awal Input -> Hidden
# ----------------------------------------------------------------------------

tabel_bobot_input_hidden <- data.frame(
  Input = c("Bias", colnames(X_train)),
  round(
    rbind(bobot_awal$v0, bobot_awal$v),
    3
  )
)

colnames(tabel_bobot_input_hidden)[-1] <-
  paste0("Hidden_", seq_len(n_hidden_terbaik))

cat("\n=== BOBOT AWAL INPUT -> HIDDEN ===\n")
print(tabel_bobot_input_hidden, row.names = FALSE)

write.csv(
  tabel_bobot_input_hidden,
  "Bobot_Awal_Input_Hidden_NNAR.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 7.2 Bobot Awal Hidden -> Output
# ----------------------------------------------------------------------------

tabel_bobot_hidden_output <- data.frame(
  Hidden = c(
    "Bias",
    paste0("Hidden_", seq_len(n_hidden_terbaik))
  ),
  Bobot = round(
    c(bobot_awal$w0, bobot_awal$w),
    3
  )
)

cat("\n=== BOBOT AWAL HIDDEN -> OUTPUT ===\n")
print(tabel_bobot_hidden_output, row.names = FALSE)

write.csv(
  tabel_bobot_hidden_output,
  "Bobot_Awal_Hidden_Output_NNAR.csv",
  row.names = FALSE
)


# ============================================================================
# 8. PELATIHAN NNAR DENGAN BOBOT AWAL
#    ARSITEKTUR 4-2-1
# ============================================================================

# ----------------------------------------------------------------------------
# 8.1 Bobot Awal
# ----------------------------------------------------------------------------

v0 <- c(0.051, -0.043)

v <- matrix(
  c(
    -0.212,  0.440,
    0.288, -0.454,
    -0.091,  0.028,
    0.383,  0.392
  ),
  nrow = 4,
  byrow = TRUE
)

w0 <- 0.178

w <- matrix(
  c(0.457, -0.047),
  nrow = 2,
  ncol = 1
)


# ----------------------------------------------------------------------------
# 8.2 Pelatihan Backpropagation
# ----------------------------------------------------------------------------

learning_rate <- lr_terbaik
max_epoch <- 10000

for (epoch in 1:max_epoch) {
  
  for (i in seq_len(nrow(X_train))) {
    
    x <- as.numeric(X_train[i, ])
    t <- Y_train[i]
    
    # Feed forward
    z <- sigmoid(v0 + x %*% v)
    y <- sigmoid(w0 + z %*% w)
    
    # Backpropagation
    error <- t - y
    
    delta_output <- error *
      turunan_sigmoid(as.numeric(y))
    
    delta_hidden <- as.numeric(delta_output) *
      as.numeric(w) *
      turunan_sigmoid(as.numeric(z))
    
    # Update bobot
    w <- w +
      learning_rate *
      as.numeric(delta_output) *
      as.numeric(z)
    
    w0 <- w0 +
      learning_rate *
      as.numeric(delta_output)
    
    v <- v +
      learning_rate *
      outer(x, delta_hidden)
    
    v0 <- v0 +
      learning_rate *
      delta_hidden
  }
}


# ============================================================================
# 9. PREDIKSI DAN EVALUASI NNAR
# ============================================================================

prediksi_norm <- numeric(nrow(X_train))

for (i in seq_len(nrow(X_train))) {
  
  x <- as.numeric(X_train[i, ])
  
  z <- sigmoid(
    v0 + x %*% v
  )
  
  prediksi_norm[i] <- sigmoid(
    w0 + z %*% w
  )
}

# Denormalisasi
prediksi_asli <- denormalisasi_minmax(
  prediksi_norm,
  y_min_train,
  y_max_train
)

aktual_asli <- denormalisasi_minmax(
  Y_train,
  y_min_train,
  y_max_train
)

# Evaluasi
rmse_nnar <- hitung_rmse(
  aktual_asli,
  prediksi_asli
)

mape_nnar <- hitung_mape(
  aktual_asli,
  prediksi_asli
)

cat("\n=== HASIL PELATIHAN NNAR ===\n")
cat("Arsitektur    : 4-2-1\n")
cat("Learning rate :", learning_rate, "\n")
cat("Epoch         :", max_epoch, "\n")
cat("RMSE          :", round(rmse_nnar, 4), "\n")
cat("MAPE          :", round(mape_nnar, 4), "%\n")


# ============================================================================
# 10. BOBOT AKHIR NNAR
# ============================================================================

# ----------------------------------------------------------------------------
# 10.1 Bobot Akhir Input -> Hidden
# ----------------------------------------------------------------------------

tabel_bobot_akhir_input_hidden <- data.frame(
  Input = c("Bias", colnames(X_train)),
  round(
    rbind(v0, v),
    3
  )
)

colnames(tabel_bobot_akhir_input_hidden)[-1] <-
  c("Hidden_1", "Hidden_2")

cat("\n=== BOBOT AKHIR INPUT -> HIDDEN ===\n")
print(
  tabel_bobot_akhir_input_hidden,
  row.names = FALSE
)

write.csv(
  tabel_bobot_akhir_input_hidden,
  "Bobot_Akhir_Input_Hidden_NNAR.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 10.2 Bobot Akhir Hidden -> Output
# ----------------------------------------------------------------------------

tabel_bobot_akhir_hidden_output <- data.frame(
  Hidden = c(
    "Bias",
    "Hidden_1",
    "Hidden_2"
  ),
  Bobot = round(
    c(w0, w),
    3
  )
)

cat("\n=== BOBOT AKHIR HIDDEN -> OUTPUT ===\n")
print(
  tabel_bobot_akhir_hidden_output,
  row.names = FALSE
)

write.csv(
  tabel_bobot_akhir_hidden_output,
  "Bobot_Akhir_Hidden_Output_NNAR.csv",
  row.names = FALSE
)


# ============================================================================
# 11. VISUALISASI ARSITEKTUR NNAR 4-2-1
# ============================================================================

# Bobot optimal
v0 <- c(0.554, 2.226)

v <- matrix(
  c(
    4.098, -4.188,
    1.800,  0.252,
    0.573,  1.240,
    0.429,  0.012
  ),
  nrow = 4,
  byrow = TRUE
)

w0 <- 0.556
w <- c(4.498, -6.973)

# Posisi neuron
yin <- c(0.80, 0.60, 0.40, 0.20)
yh <- c(0.65, 0.35)
yo <- 0.50

# Pengaturan tampilan
warna_neuron <- "#F8F4E8"
warna_garis  <- "#7A7A7A"
warna_bias   <- "#B8860B"
warna_judul  <- "#17365D"

par(mar = c(1, 1, 3, 1))

plot(
  NA,
  xlim = c(0.3, 3.7),
  ylim = c(0, 1),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  bty = "n",
  main = "Arsitektur Jaringan NNAR (4-2-1)",
  col.main = warna_judul,
  font.main = 2,
  cex.main = 1.3
)

# Input -> Hidden
for (i in 1:4) {
  for (j in 1:2) {
    
    segments(
      1, yin[i], 2, yh[j],
      col = warna_garis,
      lwd = 1.2
    )
    
    text(
      1.5,
      (yin[i] + yh[j]) / 2,
      format(v[i, j], nsmall = 3),
      cex = 0.68,
      col = warna_judul
    )
  }
}

# Hidden -> Output
for (j in 1:2) {
  
  segments(
    2, yh[j], 3, yo,
    col = warna_garis,
    lwd = 1.5
  )
  
  text(
    2.5,
    (yh[j] + yo) / 2,
    format(w[j], nsmall = 3),
    cex = 0.70,
    col = warna_judul
  )
}

# Neuron input
points(
  rep(1, 4), yin,
  pch = 21,
  bg = warna_neuron,
  col = warna_judul,
  lwd = 2,
  cex = 2.7
)

text(
  0.88, yin,
  c(
    expression(Y[t-1]),
    expression(Y[t-2]),
    expression(Y[t-3]),
    expression(Y[t-12])
  ),
  pos = 2,
  cex = 0.85,
  col = warna_judul
)

# Neuron hidden
points(
  rep(2, 2), yh,
  pch = 21,
  bg = warna_neuron,
  col = warna_bias,
  lwd = 2.5,
  cex = 2.9
)

text(
  2, yh,
  c("H1", "H2"),
  cex = 0.75,
  font = 2,
  col = warna_judul
)

# Neuron output
points(
  3, yo,
  pch = 21,
  bg = warna_neuron,
  col = warna_bias,
  lwd = 2.5,
  cex = 3
)

text(
  3, yo,
  expression(hat(Y)[t]),
  cex = 0.85,
  font = 2,
  col = warna_judul
)

# Bias hidden
points(
  1.55, 0.94,
  pch = 21,
  bg = warna_neuron,
  col = warna_bias,
  lwd = 2,
  cex = 2
)

text(
  1.55, 0.94,
  "B",
  font = 2,
  cex = 0.7,
  col = warna_bias
)

for (j in 1:2) {
  
  segments(
    1.55, 0.94,
    2, yh[j],
    col = warna_bias,
    lty = 2,
    lwd = 1.2
  )
  
  text(
    1.72,
    (0.94 + yh[j]) / 2,
    format(v0[j], nsmall = 3),
    cex = 0.65,
    col = warna_bias
  )
}

# Bias output
points(
  2.55, 0.94,
  pch = 21,
  bg = warna_neuron,
  col = warna_bias,
  lwd = 2,
  cex = 2
)

text(
  2.55, 0.94,
  "B",
  font = 2,
  cex = 0.7,
  col = warna_bias
)

segments(
  2.55, 0.94,
  3, yo,
  col = warna_bias,
  lty = 2,
  lwd = 1.2
)

text(
  2.68, 0.73,
  format(w0, nsmall = 3),
  cex = 0.65,
  col = warna_bias
)

# Label layer
text(
  1, 0.06,
  "INPUT LAYER",
  font = 2,
  cex = 0.85,
  col = warna_judul
)

text(
  2, 0.06,
  "HIDDEN LAYER",
  font = 2,
  cex = 0.85,
  col = warna_judul
)

text(
  3, 0.06,
  "OUTPUT LAYER",
  font = 2,
  cex = 0.85,
  col = warna_judul
)

# ============================================================================
# 12. PERAMALAN DAN EVALUASI NNAR(4-2-1) DATA TESTING
# ============================================================================

# ----------------------------------------------------------------------------
# 12.1 Persiapan Data Testing
# ----------------------------------------------------------------------------

y_all_norm <- c(
  y_train_norm,
  y_test_norm
)

n_train <- length(y_train_norm)
n_test  <- length(y_test_norm)

forecast_norm <- numeric(n_test)


# ----------------------------------------------------------------------------
# 12.2 Peramalan Secara Berurutan
# ----------------------------------------------------------------------------

for (h in 1:n_test) {
  
  idx <- n_train + h
  
  # Input NNAR(4-2-1)
  x <- c(
    y_all_norm[idx - 1],   # Yt-1
    y_all_norm[idx - 2],   # Yt-2
    y_all_norm[idx - 3],   # Yt-3
    y_all_norm[idx - 12]   # Yt-12
  )
  
  # Feed forward
  z <- sigmoid(
    v0 + x %*% v
  )
  
  y_hat <- sigmoid(
    w0 + z %*% w
  )
  
  forecast_norm[h] <- as.numeric(y_hat)
  
  # Hasil ramalan digunakan untuk periode berikutnya
  y_all_norm[idx] <- forecast_norm[h]
}


# ----------------------------------------------------------------------------
# 12.3 Denormalisasi Hasil Peramalan
# ----------------------------------------------------------------------------

forecast_asli <- denormalisasi_minmax(
  forecast_norm,
  y_min_train,
  y_max_train
)

aktual_testing <- y_test


# ----------------------------------------------------------------------------
# 12.4 Evaluasi
# ----------------------------------------------------------------------------

rmse_testing_nnar <- hitung_rmse(
  aktual_testing,
  forecast_asli
)

mape_testing_nnar <- hitung_mape(
  aktual_testing,
  forecast_asli
)

cat("\n=== EVALUASI NNAR(4-2-1) DATA TESTING ===\n")
cat("RMSE :", round(rmse_testing_nnar, 4), "\n")
cat("MAPE :", round(mape_testing_nnar, 4), "%\n")


# ----------------------------------------------------------------------------
# 12.5 Tabel Peramalan Data Testing
# ----------------------------------------------------------------------------

tabel_peramalan_nnar <- data.frame(
  Periode = time(y_test),
  Peramalan = round(forecast_asli, 4),
  Aktual = round(aktual_testing, 4)
)

cat("\n=== HASIL PERAMALAN NNAR(4-2-1) ===\n")
print(
  tabel_peramalan_nnar,
  row.names = FALSE
)

write.csv(
  tabel_peramalan_nnar,
  "Hasil_Peramalan_NNAR_4-2-1_Testing.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 12.6 Ringkasan Evaluasi
# ----------------------------------------------------------------------------

tabel_evaluasi_nnar <- data.frame(
  Model = "NNAR(4-2-1)",
  RMSE = round(rmse_testing_nnar, 4),
  MAPE = round(mape_testing_nnar, 4)
)

cat("\n=== RINGKASAN EVALUASI NNAR(4-2-1) ===\n")
print(
  tabel_evaluasi_nnar,
  row.names = FALSE
)

write.csv(
  tabel_evaluasi_nnar,
  "Evaluasi_NNAR_4-2-1_Testing.csv",
  row.names = FALSE
)


##############################################################################
# HYBRID SARIMA-NNAR
##############################################################################

#=============================================================================
# 1. RESIDUAL SARIMA TERBAIK
# Model: SARIMA(0,1,0)(0,1,1)[12]
#=============================================================================

theta <- sarima_terbaik$param[1]
z <- y_bc
w <- diff_manual(diff_manual(z, lag = 12), lag = 1)

residual_sarima <- numeric(length(w))
for (t in seq_along(w)) {
  residual_sarima[t] <- w[t] -
    ifelse(t > 12, theta * residual_sarima[t - 12], 0)
}

cat("\nJumlah residual SARIMA :", length(residual_sarima), "\n")


#=============================================================================
# 2. PACF RESIDUAL SARIMA
#=============================================================================

pacf_obj <- pacf(residual_sarima, plot = FALSE)
pacf_residual <- as.numeric(pacf_obj$acf)
lag_pacf <- seq_along(pacf_residual)
batas_pacf_residual <- 1.96 / sqrt(length(residual_sarima))

par(mar = c(4, 4, 3, 1))
plot(lag_pacf, pacf_residual, type = "h", lwd = 2,
     main = "PACF Residual SARIMA(0,1,0)(0,1,1)[12]",
     xlab = "Lag", ylab = "PACF",
     ylim = c(-1, 1))
abline(h = c(-batas_pacf_residual, batas_pacf_residual),
       lty = 2, lwd = 1.5)
abline(h = 0)
signifikan <- abs(pacf_residual) > batas_pacf_residual
points(lag_pacf[signifikan], pacf_residual[signifikan],
       pch = 19, cex = 1.1)
text(lag_pacf[signifikan], pacf_residual[signifikan],
     labels = lag_pacf[signifikan],
     pos = ifelse(pacf_residual[signifikan] > 0, 3, 1),
     cex = 0.8)
grid()
par(mar = c(5, 4, 4, 2) + 0.1)


#=============================================================================
# 3. NORMALISASI RESIDUAL SARIMA
#=============================================================================

normalisasi_minmax <- function(x, batas = c(0, 1)) {
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (xmax == xmin) stop("Nilai minimum dan maksimum residual sama.")
  data <- (x - xmin) / (xmax - xmin)
  data <- batas[1] + data * (batas[2] - batas[1])
  list(data = data, min = xmin, max = xmax)
}

denormalisasi_minmax <- function(x, xmin, xmax, batas = c(0, 1)) {
  xmin + ((x - batas[1]) / (batas[2] - batas[1])) * (xmax - xmin)
}

norm_res <- normalisasi_minmax(residual_sarima)
residual_norm <- norm_res$data
res_min <- norm_res$min
res_max <- norm_res$max

tabel_normalisasi_residual <- data.frame(
  Residual = round(residual_sarima, 4),
  Norm = round(residual_norm, 4)
)

cat("\n=== HASIL NORMALISASI RESIDUAL ===\n")
print(tabel_normalisasi_residual, row.names = FALSE)


#=============================================================================
# 4. PEMBENTUKAN INPUT NNAR
#    Input: Yt-1, Yt-2, Yt-3
#=============================================================================

buat_lag_nnar <- function(y, p = 3) {
  n <- length(y)
  if (n <= p) stop("Jumlah data tidak cukup untuk membentuk lag NNAR.")
  
  data <- data.frame(
    t = (p + 1):n,
    Yt = y[(p + 1):n]
  )
  
  for (lag in 1:p) {
    data[[paste0("Yt_", lag)]] <-
      y[(p + 1 - lag):(n - lag)]
  }
  data
}

lampiran_input <- buat_lag_nnar(residual_norm, p = 3)

input_nnar <- c("Yt_1", "Yt_2", "Yt_3")
X_res <- as.matrix(lampiran_input[, input_nnar])
Y_res <- as.numeric(lampiran_input$Yt)

tabel_input_nnar <- data.frame(
  t = lampiran_input$t,
  Yt_1 = round(lampiran_input$Yt_1, 4),
  Yt_2 = round(lampiran_input$Yt_2, 4),
  Yt_3 = round(lampiran_input$Yt_3, 4),
  Yt = round(lampiran_input$Yt, 4)
)

cat("\n=== DATA INPUT NNAR RESIDUAL ===\n")
print(tabel_input_nnar, row.names = FALSE)


#=============================================================================
# 5. FUNGSI NNAR
#    Arsitektur: 3-2-1
#=============================================================================

sigmoid <- function(x) {
  x <- pmax(pmin(x, 500), -500)
  1 / (1 + exp(-x))
}

hitung_rmse <- function(aktual, prediksi) {
  sqrt(mean((aktual - prediksi)^2, na.rm = TRUE))
}

hitung_mape <- function(aktual, prediksi) {
  indeks <- abs(aktual) > .Machine$double.eps
  if (!any(indeks)) return(NA_real_)
  mean(abs((aktual[indeks] - prediksi[indeks]) / aktual[indeks])) * 100
}

inisialisasi_bobot <- function(n_input, n_hidden, seed = 123) {
  set.seed(seed)
  v <- matrix(
    runif((n_input + 1) * n_hidden, -0.5, 0.5),
    nrow = n_input + 1,
    ncol = n_hidden
  )
  w <- runif(n_hidden + 1, -0.5, 0.5)
  list(v = v, w = w)
}

prediksi_nnar <- function(X, model) {
  Xb <- cbind(1, X)
  H <- sigmoid(Xb %*% model$v)
  Yhat <- cbind(1, H) %*% model$w
  as.numeric(Yhat)
}

latih_nnar <- function(X, Y, n_hidden, learning_rate,
                       threshold = 0.001, max_epoch = 10000,
                       seed = 123) {
  
  n_input <- ncol(X)
  bobot <- inisialisasi_bobot(n_input, n_hidden, seed)
  v <- bobot$v
  w <- bobot$w
  
  Xb <- cbind(1, X)
  n <- nrow(X)
  epoch_akhir <- max_epoch
  
  for (epoch in 1:max_epoch) {
    max_update <- 0
    
    for (i in 1:n) {
      x <- Xb[i, ]
      
      # Hidden layer
      net_h <- as.numeric(x %*% v)
      h <- sigmoid(net_h)
      
      # Output layer
      h_b <- c(1, h)
      y_hat <- sum(h_b * w)
      
      # Error
      error <- Y[i] - y_hat
      
      # Delta output
      delta_o <- error
      
      # Delta hidden
      delta_h <- h * (1 - h) * w[-1] * delta_o
      
      # Update bobot output
      update_w <- learning_rate * delta_o * h_b
      w_baru <- w + update_w
      
      # Update bobot input-hidden
      update_v <- learning_rate * outer(x, delta_h)
      v_baru <- v + update_v
      
      max_update <- max(
        max_update,
        max(abs(update_w)),
        max(abs(update_v))
      )
      
      w <- w_baru
      v <- v_baru
    }
    
    if (max_update < threshold) {
      epoch_akhir <- epoch
      break
    }
  }
  
  list(
    v = v,
    w = w,
    epoch = epoch_akhir
  )
}

#=============================================================================
# PARAMETER NNAR
#=============================================================================

threshold_nnar <- 0.1
max_epoch_nnar <- 10000
n_hidden_terbaik <- 2
learning_rate_list <- seq(0.01, 0.10, by = 0.01)

#=============================================================================
# 6. TRIAL AND ERROR LEARNING RATE
#=============================================================================

hasil_trial <- data.frame()

for (lr in learning_rate_list) {
  model <- latih_nnar(
    X = X_res, Y = Y_res,
    n_hidden = n_hidden_terbaik,
    learning_rate = lr,
    threshold = threshold_nnar,
    max_epoch = max_epoch_nnar,
    seed = 123
  )
  
  pred_norm <- prediksi_nnar(X_res, model)
  pred <- denormalisasi_minmax(pred_norm, res_min, res_max)
  aktual <- denormalisasi_minmax(Y_res, res_min, res_max)
  
  hasil_trial <- rbind(
    hasil_trial,
    data.frame(
      Input = length(input_nnar),
      Hidden = n_hidden_terbaik,
      Learning_Rate = lr,
      RMSE = hitung_rmse(aktual, pred),
      MAPE = hitung_mape(aktual, pred),
      Epoch = model$epoch
    )
  )
}

hasil_trial$RMSE <- round(hasil_trial$RMSE, 4)
hasil_trial$MAPE <- round(hasil_trial$MAPE, 4)

cat("\n=== HASIL TRIAL AND ERROR NNAR RESIDUAL ===\n")
print(hasil_trial, row.names = FALSE)


#=============================================================================
# 7. PEMILIHAN LEARNING RATE
#=============================================================================

# Pemilihan berdasarkan MAPE terkecil
indeks_terbaik <- which.min(hasil_trial$MAPE)
terbaik_nnar <- hasil_trial[indeks_terbaik, ]
lr_terbaik <- terbaik_nnar$Learning_Rate
n_hidden_terbaik <- terbaik_nnar$Hidden

cat("\n=== MODEL NNAR RESIDUAL TERBAIK ===\n")
cat("Input neuron  :", terbaik_nnar$Input, "\n")
cat("Hidden neuron :", terbaik_nnar$Hidden, "\n")
cat("Learning rate :", format(lr_terbaik, nsmall = 2), "\n")
cat("RMSE          :", terbaik_nnar$RMSE, "\n")
cat("MAPE          :", terbaik_nnar$MAPE, "%\n")

tabel_learning_rate <- hasil_trial[, c(
  "Learning_Rate", "RMSE", "MAPE"
)]
names(tabel_learning_rate)[1] <- "Learning Rate"

cat("\n=== TABEL HASIL LEARNING RATE ===\n")
print(tabel_learning_rate, row.names = FALSE)


#=============================================================================
# 8. BOBOT AWAL MODEL TERBAIK
#=============================================================================

bobot_awal <- inisialisasi_bobot(
  n_input = length(input_nnar),
  n_hidden = n_hidden_terbaik,
  seed = 123
)

tabel_bobot_input_hidden <- data.frame(
  Input = c("Bias", "Yt_1", "Yt_2", "Yt_3"),
  Hidden_1 = round(bobot_awal$v[, 1], 4),
  Hidden_2 = round(bobot_awal$v[, 2], 4)
)

tabel_bobot_hidden_output <- data.frame(
  Hidden = c("Bias", "Hidden_1", "Hidden_2"),
  Bobot = round(bobot_awal$w, 4)
)

cat("\n=== BOBOT AWAL INPUT KE HIDDEN ===\n")
print(tabel_bobot_input_hidden, row.names = FALSE)

cat("\n=== BOBOT AWAL HIDDEN KE OUTPUT ===\n")
print(tabel_bobot_hidden_output, row.names = FALSE)


#=============================================================================
# 9. PELATIHAN MODEL NNAR TERBAIK
#=============================================================================

model_nnar_terbaik <- latih_nnar(
  X = X_res,
  Y = Y_res,
  n_hidden = n_hidden_terbaik,
  learning_rate = lr_terbaik,
  threshold = threshold_nnar,
  max_epoch = max_epoch_nnar,
  seed = 123
)

tabel_bobot_optimal_input_hidden <- data.frame(
  Input = c("Bias", "Yt_1", "Yt_2", "Yt_3"),
  Hidden_1 = round(model_nnar_terbaik$v[, 1], 4),
  Hidden_2 = round(model_nnar_terbaik$v[, 2], 4)
)

tabel_bobot_optimal_hidden_output <- data.frame(
  Hidden = c("Bias", "Hidden_1", "Hidden_2"),
  Bobot = round(model_nnar_terbaik$w, 4)
)

cat("\n=== BOBOT OPTIMAL INPUT KE HIDDEN ===\n")
print(tabel_bobot_optimal_input_hidden, row.names = FALSE)

cat("\n=== BOBOT OPTIMAL HIDDEN KE OUTPUT ===\n")
print(tabel_bobot_optimal_hidden_output, row.names = FALSE)

cat("\nEpoch akhir :", model_nnar_terbaik$epoch, "\n")


#=============================================================================
# 10. VISUALISASI ARSITEKTUR HYBRID SARIMA-NNAR
#=============================================================================

v <- model_nnar_terbaik$v[2:4, , drop = FALSE]
w <- model_nnar_terbaik$w[2:3]
v0 <- model_nnar_terbaik$v[1, ]
w0 <- model_nnar_terbaik$w[1]

warna_neuron <- "#F8F4E8"
warna_garis <- "#7A7A7A"
warna_bias <- "#B8860B"
warna_judul <- "#17365D"

yin <- c(0.75, 0.50, 0.25)
yh <- c(0.65, 0.35)
yo <- 0.50

par(mar = c(1, 1, 3, 1))
plot(NA, xlim = c(0.3, 3.7), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     bty = "n", main = "Arsitektur Hybrid SARIMA-NNAR (3-2-1)",
     col.main = warna_judul, font.main = 2, cex.main = 1.3)

for (i in 1:3) {
  for (j in 1:2) {
    segments(1, yin[i], 2, yh[j], col = warna_garis, lwd = 1.5)
    text(1.5, (yin[i] + yh[j]) / 2, format(v[i, j], nsmall = 4),
         cex = 0.60, col = warna_judul)
  }
}

for (j in 1:2) {
  segments(2, yh[j], 3, yo, col = warna_garis, lwd = 1.5)
  text(2.5, (yh[j] + yo) / 2, format(w[j], nsmall = 4),
       cex = 0.65, col = warna_judul)
}

points(rep(1, 3), yin, pch = 21, bg = warna_neuron,
       col = warna_judul, lwd = 2, cex = 2.7)
text(0.88, yin, c(expression(Y[t-1]), expression(Y[t-2]), expression(Y[t-3])),
     pos = 2, cex = 0.85, col = warna_judul)

points(rep(2, 2), yh, pch = 21, bg = warna_neuron,
       col = warna_bias, lwd = 2.5, cex = 3)
text(rep(2, 2), yh, c("H1", "H2"), cex = 0.75,
     font = 2, col = warna_judul)

points(3, yo, pch = 21, bg = warna_neuron,
       col = warna_bias, lwd = 2.5, cex = 3)
text(3, yo, expression(hat(e)[t]), cex = 0.85,
     font = 2, col = warna_judul)

points(1.55, 0.94, pch = 21, bg = warna_neuron,
       col = warna_bias, lwd = 2, cex = 2)
text(1.55, 0.94, "B", font = 2, cex = 0.7, col = warna_bias)

for (j in 1:2)
  segments(1.55, 0.94, 2, yh[j], col = warna_bias, lty = 2, lwd = 1.2)

points(2.55, 0.94, pch = 21, bg = warna_neuron,
       col = warna_bias, lwd = 2, cex = 2)
text(2.55, 0.94, "B", font = 2, cex = 0.7, col = warna_bias)
segments(2.55, 0.94, 3, yo, col = warna_bias, lty = 2, lwd = 1.2)

text(1, 0.06, "INPUT LAYER", font = 2, cex = 0.85, col = warna_judul)
text(2, 0.06, "HIDDEN LAYER", font = 2, cex = 0.85, col = warna_judul)
text(3, 0.06, "OUTPUT LAYER", font = 2, cex = 0.85, col = warna_judul)

par(mar = c(5, 4, 4, 2) + 0.1)


#=============================================================================
# 11. PERAMALAN RESIDUAL NNAR SECARA REKURSIF
#=============================================================================

h_test <- length(y_test)
n_res_train <- length(residual_norm)

residual_norm_ext <- c(residual_norm, rep(NA_real_, h_test))
forecast_residual_norm <- numeric(h_test)

v0_fc <- model_nnar_terbaik$v[1, ]
v_fc <- model_nnar_terbaik$v[2:4, , drop = FALSE]
w0_fc <- model_nnar_terbaik$w[1]
w_fc <- model_nnar_terbaik$w[2:3]

for (i in 1:h_test) {
  idx <- n_res_train + i
  
  x <- c(
    residual_norm_ext[idx - 1],
    residual_norm_ext[idx - 2],
    residual_norm_ext[idx - 3]
  )
  
  h <- numeric(n_hidden_terbaik)
  for (j in 1:n_hidden_terbaik) {
    h[j] <- sigmoid(v0_fc[j] + sum(x * v_fc[, j]))
  }
  
  y_hat <- w0_fc + sum(h * w_fc)
  forecast_residual_norm[i] <- y_hat
  residual_norm_ext[idx] <- y_hat
}

forecast_residual <- denormalisasi_minmax(
  forecast_residual_norm, res_min, res_max
)

cat("\n=== PERAMALAN RESIDUAL NNAR ===\n")
print(round(forecast_residual, 6))


#=============================================================================
# 12. KOMBINASI HYBRID
#=============================================================================

ramalan_hybrid_bc <- ramalan_bc + forecast_residual


#=============================================================================
# 13. INVERS TRANSFORMASI BOX-COX
#=============================================================================

ramalan_hybrid <- ramalan_hybrid_bc

if (nrow(hasil_boxcox$ringkasan) > 1) {
  for (i in (nrow(hasil_boxcox$ringkasan) - 1):1) {
    lambda <- hasil_boxcox$ringkasan$Lambda[i + 1]
    shift <- hasil_boxcox$ringkasan$Shift[i + 1]
    ramalan_hybrid <- inverse_boxcox(ramalan_hybrid, lambda) - shift
  }
}


#=============================================================================
# 14. EVALUASI HYBRID SARIMA-NNAR
#=============================================================================

if (length(ramalan_hybrid) != length(y_test))
  stop("Panjang ramalan hybrid tidak sama dengan data testing.")

error_hybrid <- y_test - ramalan_hybrid
rmse_hybrid <- sqrt(mean(error_hybrid^2, na.rm = TRUE))
mape_hybrid <- hitung_mape(y_test, ramalan_hybrid)

cat("\n=== EVALUASI HYBRID SARIMA-NNAR DATA TESTING ===\n")
cat("RMSE :", round(rmse_hybrid, 4), "\n")
cat("MAPE :", round(mape_hybrid, 4), "%\n")


#=============================================================================
# 15. TABEL HASIL PERAMALAN HYBRID
#=============================================================================

tabel_peramalan_hybrid <- data.frame(
  Periode = format(tanggal_test, "%Y-%m"),
  Peramalan = round(ramalan_hybrid, 3),
  Aktual = round(y_test, 3)
)

cat("\n=== HASIL PERAMALAN HYBRID SARIMA-NNAR ===\n")
print(tabel_peramalan_hybrid, row.names = FALSE)


#=============================================================================
# 16. RINGKASAN EVALUASI
#=============================================================================

tabel_evaluasi_hybrid <- data.frame(
  Model = "Hybrid SARIMA-NNAR (3-2-1)",
  RMSE = round(rmse_hybrid, 4),
  MAPE = round(mape_hybrid, 4)
)

cat("\n=== RINGKASAN EVALUASI HYBRID SARIMA-NNAR ===\n")
print(tabel_evaluasi_hybrid, row.names = FALSE)


# ============================================================================
# PERBANDINGAN AKHIR KETIGA MODEL PADA DATA TESTING
# ============================================================================

tabel_perbandingan_akhir <- data.frame(
  Model = c("SARIMA(0,1,0)(0,1,1)[12]",
            "NNAR(4-2-1)",
            "Hybrid (SARIMA(0,1,0)(0,1,1)[12]-NNAR(3-2-1))"),
  RMSE  = round(c(RMSE, rmse_testing_nnar, rmse_hybrid), 4),
  MAPE  = round(c(MAPE, mape_testing_nnar, mape_hybrid), 4)
)

cat("\n=== PERBANDINGAN AKHIR KETIGA MODEL (DATA TESTING) ===\n")
print(tabel_perbandingan_akhir, row.names = FALSE)

write.csv(
  tabel_perbandingan_akhir,
  "Perbandingan_Akhir_Model.csv",
  row.names = FALSE
)

model_terbaik_akhir <- tabel_perbandingan_akhir$Model[
  which.min(tabel_perbandingan_akhir$RMSE)
]

cat("\nModel terbaik berdasarkan RMSE terendah pada data testing:",
    model_terbaik_akhir, "\n")

# ============================================================================
# PERAMALAN 12 BULAN KE DEPAN (JANUARI - DESEMBER 2026)
# MENGGUNAKAN MODEL TERBAIK: NNAR(4-2-1)
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Bobot Optimal NNAR(4-2-1)
#    (sesuai bobot akhir yang digunakan pada evaluasi data testing)
# ----------------------------------------------------------------------------

v0 <- c(0.554, 2.226)

v <- matrix(
  c(
    4.098, -4.188,
    1.800,  0.252,
    0.573,  1.240,
    0.429,  0.012
  ),
  nrow = 4,
  byrow = TRUE
)

w0 <- 0.556
w  <- c(4.498, -6.973)

# ----------------------------------------------------------------------------
# 2. Normalisasi Seluruh Data (Training + Testing)
#    Menggunakan parameter normalisasi dari data training (y_min_train, y_max_train)
#    agar konsisten dengan skala yang dipakai saat pelatihan model
# ----------------------------------------------------------------------------

y_full_norm <- (ntp - y_min_train) / (y_max_train - y_min_train)

n_full <- length(y_full_norm)
h_future <- 12

y_full_norm_ext <- c(y_full_norm, rep(NA, h_future))

# ----------------------------------------------------------------------------
# 3. Peramalan Rekursif 12 Bulan ke Depan
# ----------------------------------------------------------------------------

forecast_future_norm <- numeric(h_future)

for (h in 1:h_future) {
  
  idx <- n_full + h
  
  x <- c(
    y_full_norm_ext[idx - 1],   # Yt-1
    y_full_norm_ext[idx - 2],   # Yt-2
    y_full_norm_ext[idx - 3],   # Yt-3
    y_full_norm_ext[idx - 12]   # Yt-12
  )
  
  z <- sigmoid(v0 + x %*% v)
  y_hat <- sigmoid(w0 + z %*% w)
  
  forecast_future_norm[h] <- as.numeric(y_hat)
  y_full_norm_ext[idx]    <- forecast_future_norm[h]
}

# ----------------------------------------------------------------------------
# 4. Denormalisasi Hasil Peramalan
# ----------------------------------------------------------------------------

forecast_future <- denormalisasi_minmax(
  forecast_future_norm,
  y_min_train,
  y_max_train
)

# ----------------------------------------------------------------------------
# 5. Tabel Hasil Peramalan 12 Bulan ke Depan
# ----------------------------------------------------------------------------

tanggal_future <- seq(
  from = tail(tanggal, 1),
  by = "month",
  length.out = h_future + 1
)[-1]

tabel_peramalan_future <- data.frame(
  Periode   = format(tanggal_future, "%Y-%m"),
  Peramalan = round(forecast_future, 3)
)

cat("\n=== PERAMALAN NTP NASIONAL 12 BULAN KE DEPAN (NNAR 4-2-1) ===\n")
print(tabel_peramalan_future, row.names = FALSE)