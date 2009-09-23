LengthIncr <-function(gr.alpha, gr.beta, L1, delta.T, alpha = 35, beta = 55){
    # alpha must be smaller than beta
    if(alpha >= beta) stop("Error: parameter alpha must be smaller than beta")
       ((beta * gr.alpha - alpha * gr.beta)/(gr.alpha - gr.beta) - L1) * (1 - ( 1 + (gr.alpha - gr.beta) / (alpha - beta)) ^ delta.T)
}

