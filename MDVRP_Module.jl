module MDVRP_Module
using Distances, DataFrames, Plots, Colors, Random, JuMP, HiGHS, Statistics,TravelingSalesmanHeuristics

export Depot, Customer, MDVRPInstance
export parse_mdvrp, instance_to_dataframe, instance_summary
export plot_instance
export create_static_df, create_dynamic_df, plot_clustering
export create_grouping_df, plot_groups
export prepare_routing_data
export solve_tsp_nnh, solve_tsp_milp, solve_milp
export compare_nnh_milp
export overlay_routes!

# ── Structs ────────────────────────────────────────────────────────────────────
struct Depot
    id    :: Int        # global index (n+1 … n+t)
    x     :: Float64
    y     :: Float64
    D     :: Float64    # max route duration
    Q     :: Float64    # max vehicle capacity
end

struct Customer
    id    :: Int        # 1 … n
    x     :: Float64
    y     :: Float64
    d     :: Float64    # service duration
    q     :: Float64    # demand
end

struct MDVRPInstance
    problem_type :: Int          # 2 = MDVRP
    k_s           :: Int   # vehicles per depot (same for all depots in this instance)
    n            :: Int          # number of customers
    t            :: Int          # number of depots
    customers    :: Vector{Customer}
    depots       :: Vector{Depot}
    coords       :: Matrix{Float64}   # 2×(n+t) — all vertices, cols = [customers | depots]
    C            :: Matrix{Float64}   # (n+t)×(n+t) pairwise Euclidean distance matrix
end
# ── Functions ──────────────────────────────────────────────────────────────────
function parse_mdvrp(filepath::String) :: MDVRPInstance

    lines = filter(!isempty, readlines(filepath))

    # Line 1: type k_s n t
    header          = parse.(Int, split(lines[1]))
    problem_type, k_s, n, t = header[1], header[2], header[3], header[4]

    # Lines 2 … t+1: depot capacities (D Q) — stored temporarily, matched to
    # depot coordinates read later from the last t vertex lines
    depot_DQ = Matrix{Float64}(undef, t, 2)
    for s in 1:t
        row          = parse.(Float64, split(lines[1 + s]))
        depot_DQ[s, :] = row[1:2]   # [D, Q]
    end

    # Lines t+2 … t+1+n+t: all vertices (customers first, then depots)
    # format: i  x  y  d  q  f  a  list...
    customers = Vector{Customer}()
    depots    = Vector{Depot}()

    for k in 1:(n + t)
        row = parse.(Float64, split(lines[1 + t + k]))
        i, x, y, d, q = Int(row[1]), row[2], row[3], row[4], row[5]

        if k <= n
            push!(customers, Customer(i, x, y, d, q))
        else
            s = k - n   # depot index 1…t
            push!(depots, Depot(i, x, y, depot_DQ[s, 1], depot_DQ[s, 2]))
        end
    end

    # Build 2×(n+t) coordinate matrix — customers first, then depots
    coords = Matrix{Float64}(undef, 2, n + t)
    for c in customers
        coords[:, c.id] = [c.x, c.y]
    end
    for (s, dep) in enumerate(depots)
        coords[:, n + s] = [dep.x, dep.y]
    end

    # (n+t)×(n+t) pairwise Euclidean distance matrix
    C = pairwise(Euclidean(), coords, dims=2)

    return MDVRPInstance(problem_type, k_s, n, t, customers, depots, coords, C)
end


function instance_to_dataframe(inst::MDVRPInstance) :: DataFrame

    # Customers
    customer_rows = DataFrame(
        id   = [c.id for c in inst.customers],
        x    = [c.x  for c in inst.customers],
        y    = [c.y  for c in inst.customers],
        q    = [c.q  for c in inst.customers],
        type = fill(:customer, inst.n)
    )

    # Depots
    depot_rows = DataFrame(
        id   = [d.id for d in inst.depots],
        x    = [d.x  for d in inst.depots],
        y    = [d.y  for d in inst.depots],
        q    = fill(0.0, inst.t),          # depots have no demand
        type = fill(:depot, inst.t)
    )

    return vcat(customer_rows, depot_rows)
end


function plot_instance(df::DataFrame, filepath::String)
    customers = df[df.type .== :customer, :]
    depots    = df[df.type .== :depot,    :]

    p = scatter(customers.x, customers.y,
        zcolor     = customers.q,
        colorbar   = true,
        label      = "Customers",
        marker     = :circle,
        markersize = 6,
        #markerstrokewidth = 0.2,      # Makes the border slimmer
        #markerstrokealpha = 1.0,      # Optional: makes border slightly transparent
        legend_background_color = RGBA(1, 1, 1, 0.8),   # ← add this
        seriescolor     = cgrad(:inferno, rev = true), 
        xlabel     = "x",
        ylabel     = "y",
        title      = basename(filepath) * " — MDVRP instance"
    )

    scatter!(p, depots.x, depots.y,
        label      = "Depots",
        marker     = :star5,
        markersize = 10,
        color      = :red
    )

    return p
end


function create_static_df(inst)

    # Step 1 ── Slice: first n rows (customers), last t columns (depots)
    sub = inst.C[1:inst.n, end-inst.t+1:end]

    # Step 2 ── Store in a DataFrame
    depot_ids = [d.id for d in inst.depots]
    dist_cols = [Symbol("dist_depot_$(id)") for id in depot_ids]

    static_df = DataFrame(sub, dist_cols)

    customer_ids = [c.id for c in inst.customers]
    insertcols!(static_df, 1, :customer => customer_ids)

    # Step 3 ── Append ranking columns
    rank_cols = ["$(r)st closest" for r in 1:inst.t]

    for col in rank_cols
        static_df[!, col] = fill(0, inst.n)
    end

    for i in 1:inst.n
        dists    = [static_df[i, col] for col in dist_cols]
        order    = sortperm(dists)
        resolved = zeros(Int, inst.t)

        r = 1
        while r <= inst.t
            current_dist = dists[order[r]]
            tied_ranks   = [s for s in r:inst.t if dists[order[s]] == current_dist]
            tied_ids     = shuffle([depot_ids[order[s]] for s in tied_ranks])

            for (k, s) in enumerate(tied_ranks)
                resolved[s] = tied_ids[k]
            end

            r += length(tied_ranks)
        end

        for r in 1:inst.t
            static_df[i, rank_cols[r]] = resolved[r]
        end
    end

    return static_df
end



#function create_dynamic_df(inst, static_df, df)
function create_dynamic_df(inst, static_df, df, instance_name::String = "")

    rank_cols = ["$(r)st closest" for r in 1:inst.t]
    customer_coords = df[df.type .== :customer, [:id, :x, :y, :q]]

    dynamic_df = DataFrame(
        customer       = [c.id for c in inst.customers],
        assigned_depot = [static_df[i, rank_cols[1]] for i in 1:inst.n],
        assigned_rule  = fill("Closeness", inst.n),
        x              = customer_coords.x,
        y              = customer_coords.y,
        q              = customer_coords.q
    )

    depot_capacity = Dict(d.id => inst.k_s * d.Q for d in inst.depots)
    depot_load     = Dict(d.id => 0.0 for d in inst.depots)

    for i in 1:inst.n
        customer_demand = inst.customers[i].q

        assigned = false
        for r in 1:inst.t
            candidate_depot = static_df[i, rank_cols[r]]

            if depot_load[candidate_depot] + customer_demand <= depot_capacity[candidate_depot]
                dynamic_df[i, :assigned_depot] = candidate_depot
                dynamic_df[i, :assigned_rule]  = r == 1 ? "Closeness" : "Capacity"
                depot_load[candidate_depot]   += customer_demand
                assigned = true
                break
            end
        end

        if !assigned
            @warn "Customer $(inst.customers[i].id) with demand $customer_demand could not be assigned — all depots full."
        end
    end

    # Build depot capacities summary instead of printing
    depot_capacities_df = DataFrame(
        Instance = fill(instance_name, inst.t),
        Depot    = [d.id for d in inst.depots],
        Load     = [depot_load[d.id] for d in inst.depots],
        Capacity = [depot_capacity[d.id] for d in inst.depots]
    )

    #for d in inst.depots
    #println("  Depot $(d.id): load = $(depot_load[d.id]) / capacity = $(depot_capacity[d.id])")
    #end
    
    return dynamic_df, depot_capacities_df
end


function plot_clustering(df::DataFrame, dynamic_df::DataFrame, filepath::String)
    depots    = df[df.type .== :depot,    :]
    customers = df[df.type .== :customer, :]

    # Join customers with their assigned depot from dynamic_df
    customers = innerjoin(customers, dynamic_df[:, [:customer, :assigned_depot]],
                          on = :id => :customer)

    # One color per depot using a distinguishable palette
    depot_ids  = sort(unique(customers.assigned_depot))
    colors     = distinguishable_colors(length(depot_ids), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    color_map  = Dict(id => colors[i] for (i, id) in enumerate(depot_ids))

    p = plot(
        xlabel = "x",
        ylabel = "y",
        title  = basename(filepath) * " — Clustering",
        legend_background_color = RGBA(1, 1, 1, 0.8),   # ← add this
    )

    for id in depot_ids
        cluster = customers[customers.assigned_depot .== id, :]
        scatter!(p, cluster.x, cluster.y,
            label      = "Cluster depot $id",
            marker     = :circle,
            markersize = 8,
            markerstrokewidth = 0.4,      # Makes the border slimmer
            markeralpha       = 0.8,
            color      = color_map[id]
        )
    end

    scatter!(p, depots.x, depots.y,
        label      = "Depots",
        marker     = :star5,
        markersize = 10,
        color      = :red
    )

    return p
end


#function create_grouping_df(inst::MDVRPInstance, dynamic_df::DataFrame) :: DataFrame
function create_grouping_df(inst::MDVRPInstance, dynamic_df::DataFrame,instance_name::String = "") :: Tuple{DataFrame, DataFrame}
    

    # Step 1 — join customer coordinates into the working dataframe
    working = copy(dynamic_df)
    working[!, :assigned_vehicle] = fill(0, nrow(working))

    # Step 2 — build a depot lookup and a global vehicle offset per depot
    depot_dict   = Dict(d.id => d for d in inst.depots)
    depot_offset = Dict(inst.depots[s].id => (s - 1) * inst.k_s 
                        for s in 1:inst.t)

    # Accumulator for the vehicles capacity summary
    vehicle_records = NamedTuple[]
    
    # Step 3 — angular sweep per depot
    for depot_id in unique(working.assigned_depot)
        depot       = depot_dict[depot_id]
        row_indices = findall(working.assigned_depot .== depot_id)

        xs = working.x[row_indices]
        ys = working.y[row_indices]
        qs = working.q[row_indices]

        # polar angle relative to this depot, sweep order
        angles = atan.(ys .- depot.y, xs .- depot.x)
        order  = sortperm(angles)

        vehicle_loads = zeros(inst.k_s)
        vehicle_idx   = 1

        for pos in order
            actual_row = row_indices[pos]
            demand     = qs[pos]

            # spillover: advance to next vehicle if current one is full
            while vehicle_idx <= inst.k_s && vehicle_loads[vehicle_idx] + demand > depot.Q
                vehicle_idx += 1
            end

            if vehicle_idx > inst.k_s
                @warn "Customer $(working.customer[actual_row]) could not fit in any vehicle at depot $depot_id — overflowed past vehicle $(depot_offset[depot_id] + inst.k_s)"
                continue
            end

            working[actual_row, :assigned_vehicle] = depot_offset[depot_id] + vehicle_idx
            vehicle_loads[vehicle_idx] += demand # ← main sweep closes here
        end

        # after the main sweep loop, still inside the depot block
        # second pass starts here, same indentation as the for loop above
        unassigned = [row_indices[pos] for pos in order 
              if working.assigned_vehicle[row_indices[pos]] == 0]

        # add this debug line right after building unassigned
        #println("Unassigned rows: ", unassigned)
        #println("Unassigned customers: ", working.customer[unassigned])

        for actual_row in unassigned
            demand = working.q[actual_row]
            cust_id = working.customer[actual_row]
            placed = false
            for v in 1:inst.k_s
                if vehicle_loads[v] + demand <= depot.Q
                    # find the true row index by customer id
                    true_row = findfirst(working.customer .== cust_id)
                    #working[actual_row, :assigned_vehicle] = depot_offset[depot_id] + v
                    #working[!, :assigned_vehicle][actual_row] = depot_offset[depot_id] + v
                    working[true_row, :assigned_vehicle] = depot_offset[depot_id] + v
                    vehicle_loads[v] += demand
                    placed = true
                    break
                end
            end
            if !placed
                @warn "Customer $(working.customer[actual_row]) genuinely cannot fit in any vehicle at depot $depot_id"
            end
        end # second pass loop closes

        # Build summary rows — one per vehicle at this depot
        for v in 1:inst.k_s
            push!(vehicle_records, (
                Instance = instance_name,
                Depot    = depot_id,
                Vehicle  = depot_offset[depot_id] + v,
                Load     = vehicle_loads[v],
                Capacity = depot.Q,
            ))
        end
        # diagnostics — inside depot loop, after second pass
        #for v in 1:inst.k_s
        #    println("  Depot $depot_id | Vehicle $(depot_offset[depot_id] + v): load = $(vehicle_loads[v]) / $(depot.Q)")
        #end
    end
    vehicles_df = DataFrame(vehicle_records)
    return working, vehicles_df
end

function plot_groups(df::DataFrame, grouping::DataFrame, filepath::String)
    depots    = df[df.type .== :depot, :]
    depot_ids = sort(unique(grouping.assigned_depot))
    colors    = distinguishable_colors(length(depot_ids), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    color_map = Dict(id => colors[i] for (i, id) in enumerate(depot_ids))

    markers = [:circle, :rect, :diamond, :utriangle, :dtriangle,
               :pentagon, :hexagon, :star4, :star6, :cross]

    p = plot(xlabel = "x", ylabel = "y",
             title  = basename(filepath) * " — Grouping",
             legend_background_color = RGBA(1, 1, 1, 0.8))

    for depot_id in depot_ids
        depot_customers = grouping[grouping.assigned_depot .== depot_id, :]
        vehicle_ids     = sort(unique(depot_customers.assigned_vehicle))

        for (i, v_id) in enumerate(vehicle_ids)
            group = depot_customers[depot_customers.assigned_vehicle .== v_id, :]

            hsv       = convert(HSV, color_map[depot_id])
            lightness = 0.6 + 0.4 * (i - 1) / max(1, length(vehicle_ids) - 1)
            color_val = HSV(hsv.h, hsv.s, lightness)

            scatter!(p, group.x, group.y,
                label             = i == 1 ? "Depot $depot_id" : "",
                color             = color_val,
                marker            = markers[mod1(i, length(markers))],
                markeralpha       = 0.8,
                markersize        = 8,
                markerstrokewidth = 0.5
            )
        end
    end

    scatter!(p, depots.x, depots.y,
        label             = "Depots",
        marker            = :star5,
        markersize        = 10,
        markerstrokewidth = 0.5,
        color             = :red
    )

    return p
end

function instance_summary(inst::MDVRPInstance)
    println("Problem type : ", inst.problem_type, " (MDVRP)")
    println("Vehicles/depot: ", inst.k_s)
    println("Customers (n): ", inst.n)
    println("Depots    (t): ", inst.t)
    println("coords size  : ", size(inst.coords))
    println("C size       : ", size(inst.C))
    println()
    println("First customer: id=$(inst.customers[1].id)  x=$(inst.customers[1].x)  y=$(inst.customers[1].y)  q=$(inst.customers[1].q)")
    println("Last depot:     id=$(inst.depots[end].id)    x=$(inst.depots[end].x)   y=$(inst.depots[end].y)   D=$(inst.depots[end].D)  Q=$(inst.depots[end].Q)")
    println()
    println("C[1,2]               = ", round(inst.C[1,2], digits=2))
    #println("C[1,161] = ", round(inst.C[1,161], digits=2), "  (distance v1 → depot 161)")
    println("C[1, first depot]    = ", round(inst.C[1, inst.n+1], digits=2))
    println("C[1, last depot]     = ", round(inst.C[1, inst.n+inst.t], digits=2))
end

function prepare_routing_data(df::DataFrame, grouping::DataFrame,
                               instance_name::String = "") :: Tuple{DataFrame, Dict{Int, NamedTuple}}

    # ── 1. Build depot coordinate lookup ──────────────────────────────────────
    depot_rows   = filter(r -> r.type == :depot, df)
    depot_coords = Dict{Int, Tuple{Float64, Float64}}(
        r.id => (r.x, r.y) for r in eachrow(depot_rows)
    )

    # ── 2. Iterate over every (depot, vehicle) group ──────────────────────────
    vehicle_groups = groupby(grouping, [:assigned_depot, :assigned_vehicle])

    routing_dict = Dict{Int, NamedTuple}()
    records      = NamedTuple[]

    for g in vehicle_groups
        depot_id   = first(g.assigned_depot)
        vehicle_id = first(g.assigned_vehicle)

        dx, dy = depot_coords[depot_id]

        cust_ids = Vector{Int}(g.customer)
        node_ids = vcat(depot_id, cust_ids)

        xs     = vcat(dx,  Vector{Float64}(g.x))
        ys     = vcat(dy,  Vector{Float64}(g.y))
        qs     = vcat(0.0, Vector{Float64}(g.q))
        coords = Matrix{Float64}(hcat(xs, ys)')
        D      = pairwise(Euclidean(), coords, dims=2)

        routing_dict[vehicle_id] = (
            depot_id   = depot_id,
            vehicle_id = vehicle_id,
            node_ids   = node_ids,
            coords     = coords,
            demands    = qs,
            D          = D,
        )

        push!(records, (
            instance    = instance_name,
            depot_id    = depot_id,
            vehicle_id  = vehicle_id,
            n_customers = length(cust_ids),
            total_load  = sum(g.q),
            node_ids    = node_ids,
            D           = D,
        ))
    end

    summary_df = DataFrame(records)
    return summary_df, routing_dict
end


function solve_tsp_nnh(routing_dict::Dict{Int, NamedTuple}) :: DataFrame

    tsp = DataFrame(
        vehicle_id = Int[],
        depot_id   = Int[],
        cost       = Float64[],
        route      = Vector{Int}[]
    )

    for (vid, v) in routing_dict
        path, cost = nearest_neighbor(v.D, firstcity=1, closepath=true)
        push!(tsp, (
            vehicle_id = vid,
            depot_id   = v.depot_id,
            cost       = round(cost, digits=4),
            route      = v.node_ids[path]
        ))
    end

    sort!(tsp, :vehicle_id)
    return tsp
end

function solve_milp(routing_dict::Dict{Int, NamedTuple}) :: DataFrame

    milp = DataFrame(
        vehicle_id = Int[],
        depot_id   = Int[],
        milp_cost  = Float64[],
        milp_route = Vector{Int}[]
    )

    for (vid, v) in routing_dict
        path, cost, status = solve_tsp_milp(v.D)

        if isnothing(path)
            @warn "Vehicle $vid — solver returned $status, skipping"
            continue
        end

        push!(milp, (
            vehicle_id = vid,
            depot_id   = v.depot_id,
            milp_cost  = round(cost, digits=4),
            milp_route = v.node_ids[path],
        ))
    end

    sort!(milp, :vehicle_id)
    return milp
end


function solve_tsp_milp(D::Matrix{Float64}) :: Tuple{Union{Vector{Int}, Nothing}, Float64, Any}
    n = size(D, 1)

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # ── Variables ─────────────────────────────────────────────────────────────
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u[2:n] >= 1)

    for i in 1:n
        fix(x[i, i], 0; force=true)
    end

    # ── Objective ─────────────────────────────────────────────────────────────
    @objective(model, Min,
        sum(D[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))

    # ── Degree constraints ────────────────────────────────────────────────────
    @constraint(model, enter[j in 1:n],
        sum(x[i, j] for i in 1:n if i != j) == 1)

    @constraint(model, leave[i in 1:n],
        sum(x[i, j] for j in 1:n if j != i) == 1)

    # ── MTZ subtour elimination ───────────────────────────────────────────────
    @constraint(model, mtz[i in 2:n, j in 2:n; i != j],
        u[i] - u[j] + n * x[i, j] <= n - 1)

    @constraint(model, [i in 2:n], u[i] <= n - 1)

    # ── Solve ─────────────────────────────────────────────────────────────────
    optimize!(model)

    status = termination_status(model)
    status != MOI.OPTIMAL && return nothing, Inf, status

    # ── Reconstruct closed tour starting and ending at depot (index 1) ────────
    x_val   = round.(Int, value.(x))
    path    = [1]
    current = 1

    for _ in 1:(n - 1)
        next = findfirst(j -> x_val[current, j] == 1, 1:n)
        isnothing(next) && break
        push!(path, next)
        current = next
        current == 1 && break
    end

    path[end] != 1 && push!(path, 1)

    return path, objective_value(model), status
end

#function compare_nnh_milp(tsp::DataFrame, milp::DataFrame) :: DataFrame
#
#    comparison = innerjoin(
#        rename(tsp[:, [:vehicle_id, :depot_id, :cost, :route]],
#               :cost  => :nnh_cost,
#               :route => :nnh_route),
#        milp[:, [:vehicle_id, :milp_cost, :milp_route]],
#        on = :vehicle_id
#    )
#
#    comparison[!, :gap_pct] = round.(
#        100.0 .* (comparison.nnh_cost .- comparison.milp_cost) ./ comparison.milp_cost,
#        digits = 2
#    )
#
#    total_nnh  = sum(comparison.nnh_cost)
#    total_milp = sum(comparison.milp_cost)
#    avg_gap    = round(mean(comparison.gap_pct), digits=2)
#
#    println("─────────────────────────────────────────────────────────────")
#    println("  NNH total cost  : ", round(total_nnh,  digits=4))
#    println("  MILP total cost : ", round(total_milp, digits=4))
#    println("  Overall gap     : ", round(100*(total_nnh - total_milp)/total_milp, digits=2), " %")
#    println("  Mean gap/vehicle: ", avg_gap, " %")
#    println("─────────────────────────────────────────────────────────────")
#    println(comparison[:, [:vehicle_id, :depot_id, :nnh_cost, :milp_cost, :gap_pct]])
#
#    return comparison
#end

function compare_nnh_milp(tsp::DataFrame, milp::DataFrame,
                           instance_name::String = "") :: Tuple{DataFrame, NamedTuple}

    comparison = innerjoin(
        rename(tsp[:, [:vehicle_id, :depot_id, :cost, :route]],
               :cost  => :nnh_cost,
               :route => :nnh_route),
        milp[:, [:vehicle_id, :milp_cost, :milp_route]],
        on = :vehicle_id
    )

    comparison[!, :gap_pct] = round.(
        100.0 .* (comparison.nnh_cost .- comparison.milp_cost) ./ comparison.milp_cost,
        digits = 2
    )

    total_nnh   = sum(comparison.nnh_cost)
    total_milp  = sum(comparison.milp_cost)
    overall_gap = round(100 * (total_nnh - total_milp) / total_milp, digits=2)
    avg_gap     = round(mean(comparison.gap_pct), digits=2)

    println("─────────────────────────────────────────────────────────────")
    println("  NNH total cost  : ", round(total_nnh,  digits=4))
    println("  MILP total cost : ", round(total_milp, digits=4))
    println("  Overall gap     : ", overall_gap, " %")
    println("  Mean gap/vehicle: ", avg_gap, " %")
    println("─────────────────────────────────────────────────────────────")
    println(comparison[:, [:vehicle_id, :depot_id, :nnh_cost, :milp_cost, :gap_pct]])

    summary_row = (
    Instance        = instance_name,
    NNH_total_cost  = total_nnh,
    MILP_total_cost = total_milp,
    Overall_gap     = overall_gap,
    Mean_gap        = avg_gap
    )
    return comparison, summary_row
end

function overlay_routes!(p, routes::DataFrame, df::DataFrame)
    # Fast id → (x,y) lookup
    coord = Dict(row.id => (row.x, row.y) for row in eachrow(df))

    for row in eachrow(routes)
        ids = row.route                        # e.g. [101, 5, 85, 91, 44, 101]
        xs  = [coord[id][1] for id in ids]
        ys  = [coord[id][2] for id in ids]

        plot!(p, xs, ys,
            label            = "",             # suppress per-segment legend clutter
            linewidth        = 1.5,
            linealpha        = 0.7,
            color            = :black,
            linestyle        = :solid,
            legend_background_color = RGBA(1, 1, 1, 0.4)
        )
    end
    return p
    end
end # module