using DataFrames

# This function aggregates df over columns given by `agg_col`. 
# The weight of each row is given in `counter_col`. 
# Columns that need to be added are given in `to_add`. 
# Columns to be averaged with or without calculating std are given in `to_mean_with_std` or `to_mean_without_std`.
# Columns not to be touched are given in `not_touch`.
# Function tries to obtain some of the above itself. 
function aggregate_data(df::DataFrame; 
    agg_col = ["G", "attr_name", "attr_degeneracy", "gamma"], counter_col = "zmax", to_add = [], 
    to_mean_with_std = [], to_mean_without_std = [], not_touch = [])

    all_known_cols = [agg_col..., counter_col, to_add..., to_mean_with_std..., to_mean_without_std..., not_touch...]

    for name in names(df)
        if name in all_known_cols
            continue
        elseif typeof(df[1, name]) == String
            append!(not_touch, [name])
        elseif endswith(name, "_std")
            append!(not_touch, [name])
        elseif length(unique(destab_res[:, name])) == 1
            append!(not_touch, [name])
        elseif name * "_std" in names(destab_res)
            append!(to_mean_with_std, [name])
        else
            append!(to_mean_without_std, [name])
        end
    end
    append!(to_add, [counter_col])

    #creating empty daraframe
    df_agg = deepcopy(df)
    delete!(df_agg, 1:nrow(df))

    df_imp = df[!, agg_col]
    df_imp_rows = eachrow(df_imp)

    for (i, row) in enumerate(eachrow(unique(df_imp)))
        inds = findall(x -> row == x, df_imp_rows)

        push!(df_agg, df[inds[1], :])

        if length(inds) > 1
            for col in to_add
                df_agg[i, col] = sum(df[inds, col])
            end

            counts = df[inds, counter_col]
            for col in to_mean_without_std
                means = df[inds, col]

                new_mean, _ = mean_std_Ksamples(counts, means, 0)

                df_agg[i, col] = new_mean
            end
            for col in to_mean_with_std
                means = df[inds, col]
                stds = df[inds, col*"_std"]

                new_mean, new_std = mean_std_Ksamples(counts, means, stds)

                df_agg[i, col] = new_mean
                df_agg[i, col*"_std"] = new_std
            end
        end
    end

    return df_agg
end



"""
This function aggregates all datapoints given in (ys, ys_err) that have the same value in xs. 
`inds` are indices of datapoints in dataset `data`. 
"""
function aggregate_data(data, inds, xs, ys, ys_err)
    partdata = data[inds, :]

    uni_xs = unique(xs)

    if length(xs) == length(uni_xs)
        return xs, ys, ys_err
    end

    uni_ys = zeros(length(uni_xs))
    uni_ys_err = zeros(length(uni_xs))

    for (i, x) in enumerate(uni_xs)
        xinds = findall(xs .== x)
        line_data = partdata[xinds, :]

        line_y = ys[xinds]
        line_y_err = []
        if !isempty(ys_err)
            line_y_err = ys_err[xinds]
        end

        uni_ys[i], uni_ys_err[i] =
            mean_std_Ksamples(line_data.repetitions, line_y, line_y_err)
    end

    return uni_xs, uni_ys, uni_ys_err
end

"""
Calculates std (and mean) of K samples of sizes `n`, means `m` and stds `s`. 
"""
function mean_std_Ksamples(n, m, s)
    mK = sum(n .* m) / sum(n)
    sK = sqrt((sum((n .- 1) .* s .^ 2) + sum(n .* (m .- mK) .^ 2)) / (sum(n) - 1))

    return mK, sK
end