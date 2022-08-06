`calculateOverlaps1` <- function(D, S, pD, pS, nrow, N, N_len,
                          ssq, B, overlaps, overlaps.P) {
 
 # Calculates overlaps in C:
 overlap <- NeedForSpeed1(D, S, pD, pS, nrow, N, N_len,
                          ssq, B, overlaps, overlaps.P)
 overlap

}
