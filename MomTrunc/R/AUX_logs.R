#returns v1/v2 where v2log = log2(v2) and v1 is positive
log2ratio0 = function(v1,v2log){
  2^(log2(v1) - v2log)
}
#returns v1/v2 where v2log = log2(v2)
log2ratio = function(v1,v2log){
  sign(v1)*2^(log2(abs(v1)) - v2log)
}


log2prod0 = function(v1,v2log){
  2^(log2(v1) + v2log)
}

log2prod = function(v1,v2log){
  sign(v1)*2^(log2(abs(v1)) + v2log)
}

# logEprod = function(v1,vlog){
#   sign(v1)*2^(log(abs(v1)) + vlog)
# }
#
#
# logEratio = function(v1,vlog){
#   sign(v1)*2^(log(abs(v1)) - vlog)
# }
