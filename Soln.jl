# include -
using PyPlot

# constants -
R = 8.314 # J/mol-K

function _data()

    # initialize -
    data_array = zeros(2,2)

    # fill -
    data_array[1,1] = 250   # K
    data_array[2,1] = 300   # K
    data_array[1,2] = 2.04  # kPa
    data_array[2,2] = 32.48 # kPa

    # return -
    return data_array
end

function _compute_slope_intercept(data_array)

    T1 = data_array[1,1]
    T2 = data_array[2,1]
    P1SAT = data_array[1,2]
    P2SAT = data_array[2,2]

    # compute b vector -
    bV = zeros(2)
    bV[1] = log(P1SAT)
    bV[2] = log(P2SAT)

    # compute X array -
    XA = zeros(2,2)
    XA[1,1] = 1
    XA[2,1] = 1
    XA[1,2] = (1/T1)
    XA[2,2] = (1/T2)

    # compute intercept and slope -
    theta_v = inv(XA)*bV

    # return -
    return theta_v
end

function _compute(T_K,theta_v)

    term = theta_v[1] + theta_v[2]*(1/T_K)
    return exp(term)
end

function _plot(updated_data_array,theta_v)

    xV = zeros(3)
    xV[1] = (1/updated_data_array[1,1])
    xV[2] = (1/updated_data_array[2,1])
    xV[3] = (1/updated_data_array[3,1])

    yV = zeros(3)
    yV[1] = log(updated_data_array[1,2])
    yV[2] = log(updated_data_array[2,2])
    yV[3] = log(updated_data_array[3,2])

    # plot data -
    plot(xV,yV,"o")

    # plot model -
    T_array = collect(200:0.1:320)
    TI_array = Float64[]
    P_array = Float64[]
    for T_value in T_array
        P_value = log(_compute(T_value,theta_v))
        push!(P_array,P_value)
        push!(TI_array,1.0/T_value)
    end

    plot(TI_array,P_array,"r--")
    xlabel(L"$\frac{1}{T}$ (K$^{-1}$)",fontsize=16)
    ylabel(L"$\ln~P_{i}^{sat}$",fontsize=16)


end

function _check(data_array)

    T1 = data_array[1,1]
    T2 = data_array[2,1]
    P1SAT = data_array[1,2]
    P2SAT = data_array[2,2]

    a = log(P2SAT/P1SAT)
    b = (1/T2 - 1/T1)
    m = (a/b)

    # compute C -
    C1 = log(P1SAT) - m*(1/T1)
    C2 = log(P2SAT) - m*(1/T2)

    return (m,C1,C2)
end

function solve()

    # load data -
    data_array = _data()

    # estimate y-intercept and slope -
    theta_v = _compute_slope_intercept(data_array)

    # compute at new T -
    T_new = 310 # K
    P_new_sat = _compute(T_new,theta_v)

    # updated data array -
    updated_data_array = [data_array ; T_new P_new_sat]

    # plot -
    _plot(updated_data_array,theta_v)

    # return -
    return (theta_v,updated_data_array)
end
