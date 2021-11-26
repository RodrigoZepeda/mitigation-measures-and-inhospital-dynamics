function fit_model(dy, y, p, t)

    #=
    Parameter specification
    =#
    beta_0   = p[1]
    beta_1   = p[2]
    beta_2   = p[3]
    beta_3   = p[4]
    beta_12  = p[5]
    beta_23  = p[6]
    eta_E    = p[7]
    gamma_0R = p[8]
    gamma_1R = p[9]
    gamma_2R = p[10]
    gamma_3R = p[11]
    theta_3M = p[12]
    p_A      = p[13]
    p_I2     = p[14]
    p_I3     = p[15]
    p_M      = p[16]
    κ1       = p[17]
    κ2       = p[18]
    κ3       = p[19]
    κ4       = p[20]
    κ5       = p[21]
    #κ6       = p[22]
    #κ7       = p[23]
    #d        = p[24]
    #t0       = p[25]

    #=
    ODE model
    =#
    #Create a better polynomial
    #=
    alpha = 1 + κ1*cos( pi*(t + t0) / d )^2 + κ2*cos( 2*pi*(t + t0) / d )^2 +
                κ3*(t/200) + κ4*(t/200)^2 + κ5*(t/200)^3 + κ6*(t/200)^4 + κ7*(t/200)^5 
    =#
    alpha = 1 + κ1*(t/200) + κ2*(t/200)^2 + κ3*(t/200)^3 + κ4*(t/200)^4 + κ5*(t/200)^5 

    #dS = -alpha(t)*(beta_0*Asymptomatic + beta_1*Infected + beta_2*UCI + 
    #         beta_3*Hospitalized)*S 
    dy[1] = -alpha*(beta_0*y[3] + beta_1*y[4] + beta_2*y[5] + beta_3*y[6])*y[1];
    
    #dE = alpha(t)*(beta_0*Asymptomatic + beta_1*Infected + beta_2*UCI + 
    #         beta_3*Hospitalized)*S - eta_E*E
    dy[2] = alpha*(beta_0*y[3] + beta_1*y[4] + beta_2*y[5] + beta_3*y[6])*y[1] - eta_E*y[2];
    
    #dA = p_A*eta_E*E - gamma_0R*A
    dy[3] = p_A*eta_E*y[2] - gamma_0R*y[3];
    
    #dI1 = (1 - p_A)*eta_E*E - (p_I2*beta_12 + (1 - p_I2)*gamma_1R)*I1
    dy[4] = (1 - p_A)*eta_E*y[2] - (p_I2*beta_12 + (1 - p_I2)*gamma_1R)*y[4];
    
    #dI2 = p_I2*beta_12*I_1 - (p_I3*beta_23 + (1 - p_I3)*gamma_2R)*I2
    dy[5] = p_I2*beta_12*y[4] - (p_I3*beta_23 + (1 - p_I3)*gamma_2R)*y[5];
    
    #dI3 = p_I3*beta_23*I2 - ((1 - pM)*gamma_3R + pM*theta_3M)*I3
    dy[6] = p_I3*beta_23*y[5] - ((1 - p_M)*gamma_3R + p_M*theta_3M)*y[6];
    
    #dM  = pM*theta_3M*I3
    dy[7] = p_M*theta_3M*y[6];

    #Cummulative infected fit
    dy[8] = (1 - p_A)*eta_E*y[2]

    #Cummulative hospitalized fit
    dy[9] = p_I2*beta_12*y[4]

    #Cummulative ICU fit
    dy[10] = p_I3*beta_23*y[5]

end
