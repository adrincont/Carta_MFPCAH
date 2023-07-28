# [Funcion para graficar objetos hibridos] ===================================.#
# [Grafica para datos funcionales multivariados] --------------------------.####
plot_func = function(FunData, name_x = 'x', name_y = 'y'){
  # FunData: funData class (funData package)
  # name_y: nombre de eje y
  # return(): ggplot class plot
  plt = list(); fig = list()
  for (i in 1:length(FunData)) {
    X = FunData[[i]]@X
    colnames(X) = FunData[[i]]@argvals[[1]]
    rownames(X) = seq(1:nrow(X))
    plt[[i]] = as.data.frame.table(X)
    colnames(plt[[i]]) = c('id', 'x', 'y')
    plt[[i]]$x = as.numeric(as.character(plt[[i]]$x))
    fig[[i]] = ggplot(plt[[i]], aes(x, y, fill = 'id')) +
      geom_line(aes(color = id), show.legend = FALSE) + theme_bw() +
      xlab(name_x) + ylab(paste(name_y, '(V', i, ')'))
  }
  if (length(fig) > 1) {return(Reduce(`+`, fig))}
  else {Reduce(`+`, fig)}
}
# Matriz de graficos para variables vectoriales
ggpairs_2 = function(data_j){
  # data: data.frame class
  # return(): ggplot class plot
  my_fn = function(data_j, mapping, ...){
    p = ggplot(data = data_j, mapping = mapping) +
      geom_point() +
      geom_smooth(method = loess, fill = 'red', color = 'red', ...) +
      geom_smooth(method = lm, fill = 'blue', color = 'blue', ...)
    p
  }
  return(ggpairs(data_j, lower = list(continuous = my_fn)))
}
# Grafica para datos funcionales hibridos
autoplot.multiFunDataH = function(object, type = c(1,2)){
  # object: multiFunDataH class
  # type: plot_func() if 1 and ggpairs_2() if 2
  # return(): ggplot class plot
  fun_data = object@funData
  vec_data = object@vecData
  # gráfico parte funcional
  list_plot = list()
  if (1 %in% type) {
    list_plot[['fun']] = plot_func(fun_data, name_x = 'value', name_y = 'y')
  }
  if (2 %in% type) {
    list_plot[['vec']] = ggpairs_2(as.data.frame(vec_data))
  }
  return(list_plot)
}
# [Fin del codigo] ===========================================================.#