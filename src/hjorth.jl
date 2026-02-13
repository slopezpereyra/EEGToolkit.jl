
"""
    hjorth_parameters(x::AbstractVector)

Calculate Hjorth's Activity, Mobility, and Complexity for a 1D signal.
Ref: Hjorth, B. (1970). EEG analysis based on time domain properties.
"""
function hjorth_parameters(x::AbstractVector)
    # 1. Activity: Variance of the signal
    activity = var(x)
    
    # First derivative approximation
    dx = diff(x)
    var_dx = var(dx)
    
    # 2. Mobility
    mobility = sqrt(var_dx / activity)
    
    # Second derivative approximation
    ddx = diff(dx)
    var_ddx = var(ddx)
    
    # 3. Complexity
    # Defined as Mobility(dx) / Mobility(x)
    mobility_dx = sqrt(var_ddx / var_dx)
    complexity = mobility_dx / mobility
    
    return (activity=activity, mobility=mobility, complexity=complexity)
end
