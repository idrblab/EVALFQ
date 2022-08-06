`calculateOverlaps2` <- function(D, pD, nrow, N, N_len,
                      B, overlaps, overlaps.P) {
 
 # Calculates overlaps in C:
 overlap <- NeedForSpeed2(D, pD, nrow, N, N_len,
                      B, overlaps, overlaps.P)
					  
 overlap
}
