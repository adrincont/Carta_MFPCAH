# [Codigos para PCA híbrido] =============================================.#####
# {Calculo para PCA funcional multiple} -----------------------------------.####
MFPCA_aux = function(mFData, M, uniExpansions, weights = rep(1, length(mFData)),
                     fit = FALSE, approx.eigen = FALSE){
  p = length(mFData)
  N = nObs(mFData)
  dimSupp = dimSupp(mFData)
  type = vapply(uniExpansions, function(l) {l$type}, FUN.VALUE = "")
  m = meanFunction(mFData, na.rm = TRUE)*0
  uniBasis = mapply(function(expansion, data) {
    do.call(univDecomp, c(list(funDataObject = data), expansion))
  }, expansion = uniExpansions, data = mFData, SIMPLIFY = FALSE)
  for (j in seq_len(p)) {
    if (type[j] == "uFPCA") 
      m[[j]] = uniBasis[[j]]$meanFunction*0
  }
  npc = vapply(uniBasis, function(x) {dim(x$scores)[2]}, FUN.VALUE = 0)
  if (M > sum(npc)) {
    M = sum(npc)
    warning("Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to ", sum(npc), ".")
  }
  if (all(foreach::foreach(j = seq_len(p), .combine = "c") %do% {uniBasis[[j]]$ortho}))
    Bchol = NULL
  else {
    Bchol = Matrix::bdiag(lapply(uniBasis, function(l) {
      if (l$ortho) 
        res = Matrix::Diagonal(n = ncol(l$scores))
      else res = Matrix::chol(l$B)
      return(res)
    }))
  }
  mArgvals = funData::argvals(mFData)
  res = MFPCA:::calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, 
                  weights = weights, npc = npc, argvals = mArgvals, uniBasis = uniBasis, 
                  fit = fit, approx.eigen = approx.eigen)
  res$meanFunction = m
  names(res$functions) = names(mFData)
  namesList = lapply(mFData, names)
  if (!all(vapply(namesList, FUN = is.null, FUN.VALUE = TRUE))) {
    if (length(unique(namesList)) != 1) 
      warning("Elements have different curve names. Use names of the first element for the results.")
    row.names(res$scores) <- namesList[[1]]
  }
  class(res) = "MFPCAfit"
  return(res)
}
# {Calculo para PCA híbrido} ----------------------------------------------.####
MFHPCA = function(datos, M = NULL, delta, 
                  omega, scale = FALSE , ...){
  if (scale) {
    datos = scale_FunDataH(datos)
  }
  fun_data = datos@funData
  vec_data = datos@vecData
  p_fun = length(fun_data) # numero de variables funciones
  p_vec =  ncol(vec_data) # numero de variables vectores
  N = nrow(vec_data) # numero de observaciones
  # {MFDA de la parte funcional}
  uniExpansions_r = list()
  for (i in 1:p_fun) {
    uniExpansions_r[[i]] = list(type = 'uFPCA',npc = 10) # PCA univariados
  }
  cat('Calculating MFPCA of the multivariate functional part...','\n')
  mfda_r  = MFPCA_aux(fun_data, M = 8,uniExpansions = uniExpansions_r)
  aport_mfda_r = summary(mfda_r)
  L = c(1:8)[round(cumsum(mfda_r$values/sum(mfda_r$values)),2) >= delta][1]
  if (is.na(L)) {
    mfda_r = MFPCA(fun_data, M = 10,uniExpansions = uniExpansions_r)
    L = c(1:10)[round(cumsum(mfda_r$values/sum(mfda_r$values)),2) >= delta][1]
  }
  if (L == 1) {
    scores_mfda = matrix(mfda_r$scores[,1:L] ,ncol = 1)
  }else{
    scores_mfda = mfda_r$scores[,1:L]
  }
  # {PCA parte vectorial}
  cat('Calculating PCA of the vectorial part...','\n')
  pca_r = PCA(vec_data, scale.unit = FALSE, ncp = p_vec, graph = FALSE)
  aport_pca_r = t(pca_r$eig)
  J = c(1:p_vec)[cumsum(pca_r$eig[,1]/sum(pca_r$eig[,1])) >= delta][1]
  pca_r = PCA(vec_data, scale.unit = FALSE, ncp = J, graph = FALSE)
  scores_pca = pca_r$ind$coord
  # {Estimacion de la matriz V}
  scores_join = cbind(scores_mfda, scores_pca)
  cat('Calculating MFHPCA of the multivariate functional part...','\n')
  V = cov(scores_join)
  eigen_V = eigen(V)
  # {Calculo de pc_hibridos}
  mfdaH_r = list()
  mfdaH_r$values = eigen_V$values # valores propios
  mfdaH_r$scores = vector() # puntajes del PCA hibrido
  cm_i = eigen_V$vectors[1:L,]
  dm_i = eigen_V$vectors[-c(1:L),]
  eta_i = scores_mfda
  gamma_i = scores_pca
  # funciones propias (Psi,Theta)
  if (is.null(M)) {
    values_mi = sort(mfdaH_r$values/sum(mfdaH_r$values), decreasing = TRUE)
    M =  c(1:(L + J))[round(cumsum(values_mi/sum(values_mi)),2) >= omega][1]
  }
  # Psi
  for (m in 1:M) {
    # Psi
    fun_i = 0
    for (l in 1:L) {
      fun_i = mfda_r$functions[l]*cm_i[l,m] + fun_i
    }
    # Theta
    vec_i = 0
    for (j in 1:J) {
      vec_i = pca_r$svd$V[,j]*dm_i[j,m] + vec_i
    }
    # Direccion hibrida
    if (m == 1) {
      dir_i = mfh_data(fun_data = fun_i,vec_data = matrix(vec_i,nrow = 1))
    }else{
      aux_dir = mfh_data(fun_data = fun_i,vec_data = matrix(vec_i,nrow = 1))
      dir_i = append_Mfun(dir_i,aux_dir)
    }
  }
  mfdaH_r$dir = dir_i
  # Calculo de scores híbridos
  for (i in 1:M) {
    D = matrix(dm_i[,i], nrow = N, ncol = J, byrow = TRUE)
    C = matrix(cm_i[,i], nrow = N, ncol = L, byrow = TRUE)
    score_i = rowSums(scores_mfda*C) + rowSums(scores_pca*D)
    mfdaH_r$scores = cbind(mfdaH_r$scores, score_i)
  }
  # Resultados finales
  aport_mfdaH_r = t(cbind(mfdaH_r$values,mfdaH_r$values/sum(mfdaH_r$values),cumsum(mfdaH_r$values/sum(mfdaH_r$values))))
  mfdaH_r$values = mfdaH_r$values[1:M]
  pca_i = list(MFPCA = mfda_r, PCA = pca_r, J = J, L = L, M = M) # PCA individuales
  pcah_i = list(V = V, mfdaH = mfdaH_r, dir = dir_i)
  aports_i = list(MFPCA = aport_mfda_r, PCA = aport_pca_r, MFPCAH = aport_mfdaH_r)
  results = list(PCAs = pca_i, MFPCAH = pcah_i,Aport = aports_i,eigen_V = eigen_V)
  class(results) = 'mfhpca'
  return(results)
}
# {Resumen PCA híbrido} ---------------------------------------------------.####
summary.mfhpca = function(object){
  mfda_r = object$PCAs$MFPCA
  pca_r = object$PCAs$PCA
  mfdaH_r = object$MFPCAH$mfdaH
  L = object$PCAs$L
  J = object$PCAs$J
  M = object$PCAs$M
  # MFPCA
  tab_summary_mfda = round(object$Aport$MFPCA, 3)[,1:L]
  tab_summary_mfda[2:3,] = tab_summary_mfda[2:3,]*100
  rownames(tab_summary_mfda) = c('val', 'explicada', 'acumulada')
  # PCA
  tab_summary_pca = round(object$Aport$PCA, 3)[,1:J]
  rownames(tab_summary_pca) = rownames(tab_summary_mfda)
  colnames(tab_summary_pca) = paste('PC', 1:J)
  # MFHPCA
  tab_summary_mfdaH = round(object$Aport$MFPCAH, 3)[,1:M]
  rownames(tab_summary_mfdaH) = rownames(tab_summary_mfda)
  colnames(tab_summary_mfdaH) = paste('PC',1:M)
  cat('|=================================================================|', '\n')
  cat('PCA _____________________________________________________________.','\n')
  print(round(tab_summary_pca,4))
  cat('MFPCA ___________________________________________________________.','\n')
  print(round(tab_summary_mfda,4))
  cat('MFHPCA __________________________________________________________.','\n')
  print(round(tab_summary_mfdaH,4))
  cat('|================================================================|','\n')
  class_scores =  class(object$MFPCAH$mfdaH$scores)
  class_dir = class(object$MFPCAH$dir)
  cat('mfdaH_r$scores:', class_scores,'...', 'puntajes para PCA hibrido','\n')
  cat('mfdaH_r$dir:', class_dir,'...',
      'direcciones principales PCA hibrido','\n')
}
# {Grafico PCA híbrido} ---------------------------------------------------.####
plot_varexp = function(values){
  val_i = values
  text_lab = paste0(round(val_i*100,1), '%')
  perf = data.frame(i = seq(1, length(val_i)), val_i, text_lab)
  max_y = max(perf$val_i)*(1 + 0.1)
  ggp = ggplot(perf)  +
    geom_bar(aes(x = i, y = val_i), stat = 'identity', fill = '#17B978') +
    geom_line(aes(x = i, y = val_i), stat = 'identity') +
    geom_point(aes(x = i, y = val_i)) +
    labs(x = 'i' ,y = 'Valor propio') +
    geom_text(aes(x = i, y = val_i), label = text_lab, vjust = -0.4, hjust = 0) +
    theme_bw() + ylim(0,max_y)
  return(ggp)
}
autoplot.mfhpca = function(object, type = 'var'){
  mfda_r = object$PCAs$MFPCA
  pca_r = object$PCAs$PCA
  mfdaH_r = object$MFPCAH$mfdaH
  L = object$PCAs$L
  J = object$PCAs$J
  M = object$PCAs$M
  k = length(object$PCAs$MFPCA$functions)
  p = dim(object$PCAs$PCA$var$coord)[1]
  N = dim(object$PCAs$PCA$ind$coord)[1]
  # Varianza explicada
  if (type == 'var') {
    p1 = plot_varexp(c(object$Aport$MFPCA[2,1:L])) + ggtitle(paste('MFPCA (k=', k,'L=', L,')'))
    p2 = plot_varexp(c(object$Aport$MFPCAH[2,1:M])) + ggtitle(paste('MFHPCA (M=' ,M ,')'))
    p3 = plot_varexp(c(object$Aport$PCA[2,1:J]/100)) + ggtitle(paste('PCA (p=', p, 'J=', J, ')'))
    plot_i = (p1 + p3)/p2
    return(plot_i)
  }
  # Gráficos de funciones propias y vectores propios PCA
  if (type == 'eigen') {
    return(autoplot(object$MFPCAH$dir))
  }
  # Plano factorial
  if (type == 'ind') {
    scores_mfdaH_r = data.frame(mfdaH_r$scores, id = 1:N)
    colnames(scores_mfdaH_r) = c('PC1', 'PC2', 'id')
    plot_i = ggplot(scores_mfdaH_r,aes(PC1, PC2)) +
      geom_vline(xintercept = 0, color = 'blue') +
      geom_hline(yintercept = 0,  color = 'blue') + theme_bw() +
      geom_point() + ggtitle('MFPCAH') +
      geom_text(aes(label = round(id,3)), vjust = -0.4, hjust = 0,size = 3) +
      geom_density_2d(alpha = 0.5)
    return(plot_i)
  }
}
# {Proyecciones con pca hibrido} ------------------------------------------.####
pryc_pca_H = function(data_j,PCA_j){
  eigen_V = PCA_j$eigen_V
  L = PCA_j$PCAs$L
  J = PCA_j$PCAs$J
  M = PCA_j$PCAs$M
  N = nObs(data_j@funData)
  cm_i = eigen_V$vectors[1:L,]
  dm_i = eigen_V$vectors[-c(1:L),]
  mfda_r = PCA_j$PCAs$MFPCA
  pca_r = PCA_j$PCAs$PCA
  scores_mfda_j = pryc_pca_F(data_j = data_j@funData,PCA_j = mfda_r)[,1:L]
  scores_pca_j = pryc_pca_V(data_j = data_j@vecData,PCA_j = pca_r)[,1:J]
  scores_j = vector()
  for (i in 1:M) {
    D = matrix(dm_i[,i], nrow = N, ncol = J, byrow = TRUE)
    C = matrix(cm_i[,i], nrow = N, ncol = L, byrow = TRUE)
    score_i = rowSums(scores_mfda_j*C) + rowSums(scores_pca_j*D)
    scores_j = cbind(scores_j, score_i)
  }
  return(scores_j)
}
pryc_pca_F = function(data_j,PCA_j){
  proy_pcaF = function(obs_w,dir_w) {
    N = nObs(obs_w)
    M = nObs(dir_w)
    # Calculo de scores híbridos
    score_i1 = matrix(nrow = N,ncol = M)
    for (m in 1:M) {
      score_i1[,m] = c(sapply(list(1:N),function(i,z_j,dir_j){
        return(scalarProduct(z_j[i],dir_j))
      },z_j = obs_w,dir_j = dir_w[m]))
    }
    colnames(score_i1) = paste0("scor_",1:M)
    return(score_i1)
  }
  Psi_i = PCA_j$functions
  scores_data = proy_pcaF(obs_w = data_j,dir_w = Psi_i)
  return(scores_data)
}
pryc_pca_V = function(data_j, PCA_j) {
  dir_j = PCA_j$svd$V
  ncp = ncol(dir_j)
  coord = data_j
  coord = crossprod(t(data_j), dir_j)
  coord = coord[, 1:ncp, drop = F]
  return(coord)
}
# {Reconstruccion de datos apartir de scors y direcciones} ----------------.####
Inv_PCA_H = function(dir_j,scor_j) {
  n_fun = length(dir_j@funData)
  n_vec = ncol(dir_j@vecData)
  N = nrow(scor_j)
  simData_fun = vector("list", n_fun)
  for (j in seq_len(n_fun)) {
    aux_mat = dir_j@funData[[j]]@X
    X <- apply(aux_mat,-1, function(v) {
      scor_j %*% v
    })
    if (N == 1) 
      dim(X) <- c(1, nObsPoints(dir_j@funData[[j]]))
    simData_fun[[j]] <- funData(dir_j@funData[[j]]@argvals, X)
  }
  simData_fun = multiFunData(simData_fun)
  simData_vec = scor_j %*% dir_j@vecData
  simData_vec = simData_vec[,1:n_vec]
  if (n_vec == 1) {
    simData_vec = matrix(c(simData_vec),ncol = 1)
  }
  Z = mfh_data(fun_data = simData_fun,vec_data = simData_vec)
  return(Z)
}
Inv_PCA_F = function(dir_j,scor_j) {
  n_fun = length(dir_j)
  N = nrow(scor_j)
  simData_fun = vector("list", n_fun)
  for (j in seq_len(n_fun)) {
    aux_mat = dir_j[[j]]@X
    X <- apply(aux_mat,-1, function(v) {
      scor_j %*% v
    })
    if (N == 1) 
      dim(X) <- c(1, nObsPoints(dir_j[[j]]))
    simData_fun[[j]] <- funData(dir_j[[j]]@argvals, X)
  }
  simData_fun = multiFunData(simData_fun)
  return(simData_fun)
}
# [Fin del codigo] ========================================================.####