# [Funcion para crear objeto híbrido] ========================================.#
# Definir clase multiFunDataH
check_FunDataH = function(object){
  errors = character()
  if (length(errors) == 0) TRUE else errors
}
setClass("multiFunDataH", representation(funData = "multiFunData", vecData = "matrix"),validity = check_FunDataH)
# Crear objeto de la clase multiFunDataH
mfh_data = function(fun_data,vec_data) {
  # fun_data: multiFunData class (funData package)
  # vec_data: matrix or array
  Xi = new("multiFunDataH", funData = fun_data,vecData = vec_data)
  return(Xi)
}
# Métodos para mostrar objeto híbrido
print_multiFunDataH = function(e1) {
  N_i = dim(e1@funData[[1]]@X)[1]
  cat('An object of class "multiFunDataH" whit N =', N_i, 'observations','\n')
  cat('_________________________________________________________________', '\n')
  cat('multiFunData:', 'k = ', length(e1@funData), '\n')
  cat('vecData:', 'p = ', ncol(e1@vecData), '\n')
  cat('_________________________________________________________________', '\n')
}
setMethod("show",signature = "multiFunDataH",function(object){print_multiFunDataH(object)})
# Método de suma entre objetos híbridos
setMethod("+", signature(e1 = "multiFunDataH", e2 = "multiFunDataH"), function(e1, e2) {
  sum_fun = e1@funData + e2@funData
  if (nrow(e1@vecData) == 1 & nrow(e2@vecData) > 1) {
    e1@vecData = matrix(e1@vecData,nrow = nrow(e2@vecData),ncol = ncol(e1@vecData),byrow = TRUE)
  } else if (nrow(e2@vecData) == 1 & nrow(e1@vecData) > 1) {
    e2@vecData = matrix(e2@vecData,nrow = nrow(e1@vecData),ncol = ncol(e2@vecData),byrow = TRUE)
  }
  sum_vec = e1@vecData + e2@vecData
  ei = mfh_data(fun_data = sum_fun,vec_data = sum_vec)
  structure(ei, class = "multiFunDataH")
})
# Metodo producto escalar para híbridos
setMethod("*", signature(e1 = "numeric", e2 = "multiFunDataH"), function(e1,e2) {
  sum_fun = e1*e2@funData
  sum_vec = e1*e2@vecData
  ei = mfh_data(fun_data = sum_fun,vec_data = sum_vec)
  structure(ei, class = "multiFunDataH")
})
# Metodo producto interno entre objetos hibridos
scalarProduct_mfh = function(Z1,Z2, w = 1){
  return(scalarProduct(Z1@funData,Z2@funData) + w*diag(Z1@vecData %*% t(Z2@vecData)))
}
# Metodo para unir observaciones para hibridos
append_Mfun = function(Z1,Z2){
  # Z1,Z2: objetos de la clase multiFunDataH que se van a unir
  # retur(): objetos de la clase multiFunDataH
  vec_append = rbind(Z1@vecData,Z2@vecData)
  n_fun = length(Z1@funData)
  append_list = list()
  fun_append = list()
  fun_append_i = rbind(Z1@funData[[1]]@X,Z2@funData[[1]]@X)
  args_append_i = Z1@funData[[1]]@argvals
  fun_append[[1]] = funData(args_append_i,fun_append_i)
  if (n_fun > 1) {
    for (j in 1:n_fun) {
      fun_append_i = rbind(Z1@funData[[j]]@X,Z2@funData[[j]]@X)
      args_append_i = Z1@funData[[j]]@argvals
      fun_append[[j]] = funData(args_append_i,fun_append_i)
    }
  }
  fmult_append = multiFunData(fun_append)
  Zi = mfh_data(fun_data = fmult_append,vec_data = vec_append)
  return(Zi)
}
# Método para seleccionar observaciones hibridas
select_Mfun = function(Z1,id){
  funData = Z1@funData[id]
  vecData = matrix(Z1@vecData[id,],ncol = dim(Z1@vecData)[2],nrow = length(id))
  return(mfh_data(fun_data = funData,vec_data = vecData))
}
# [Fin del codigo] ===========================================================.#