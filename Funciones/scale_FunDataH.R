# [Función para estandarizar objetos híbridos] ===============================.#
# Escalador min max
MinMax = function(x, args = NULL){
  if (is.null(args)) {args = c(min(x),max(x))}
  val_i = (x - args[1])/(args[2] - args[1])
  attr(val_i, 'args') = args
  return(val_i)
}
# Calcular objetos para la estandarizacion
scale_FunDataH_Train = function(Z) {
  # Z: objeto de la clase multiFunDataH
  fun_data = Z@funData
  vec_data = Z@vecData
  # Escalar variables funcionales
  scale_y = max(sapply(fun_data@.Data, function(x){max(x@X)}))
  scale_x = list()
  for (i in 1:length(fun_data@.Data)) {
    min_max_i = MinMax(fun_data@.Data[[i]]@argvals[[1]])
    scale_x[[i]] = attr(min_max_i, 'args')
  }
  scal_fun = list(scale_y = scale_y,scale_x = scale_x)
  # Media híbrida
  mean_fun = meanFunction(fun_data)
  mean_vect = apply(vec_data, 2, mean)
  mean_fh = mfh_data(fun_data = mean_fun,vec_data = matrix(mean_vect,nrow = 1))
  return(list(scal_fun = scal_fun,mean_fh = mean_fh))
}
# Estandarizacion de datos híbridos
scale_FunDataH = function(Z, scale_obj = NULL, center = TRUE, w = NULL) {
  if (is.null(scale_obj)) {
    scale_obj = scale_FunDataH_Train(Z)
  }
  # Centrar por media hibrida
  if (center) {
    Z@funData = Z@funData - scale_obj$mean_fh@funData
    mean_vec = matrix(scale_obj$mean_fh@vecData,
                      nrow = nrow(Z@vecData),ncol = ncol(Z@vecData),byrow = TRUE)
    Z@vecData = Z@vecData - mean_vec 
  }
  # Z: objeto de la clase multiFunDataH
  fun_data = Z@funData
  vec_data = Z@vecData
  # Escalar variables funcionales
  fun_data = fun_data/scale_obj$scal_fun$scale_y
  for (i in 1:length(fun_data@.Data)) {
    scale_x = scale_obj$scal_fun$scale_x[[i]]
    min_max_i = MinMax(fun_data@.Data[[i]]@argvals[[1]],args = scale_x)
    fun_data@.Data[[i]]@argvals[[1]] = min_max_i
  }
  # Eliminar efecto entre funciones y vectores
  if (center) {fun_data_c = fun_data}else{
    fun_data_c = fun_data - meanFunction(fun_data)
  }
  if (is.null(w)) {
    w = abs(sum(abs(scalarProduct(fun_data_c, fun_data_c)))/sum(rowSums(vec_data^2)))
  }
  vec_data = sqrt(w)*vec_data
  scala_fundatah = mfh_data(fun_data = fun_data,vec_data = vec_data)
  attr(scala_fundatah, 'w') = w
  return(scala_fundatah)
}
# [Fin del codigo] ===========================================================.#