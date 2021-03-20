

"""
Function: ArtificalBase_table
Argument: b: a vector of node supply and demand
Return: artificalBase_table: an inital look up table for the
    artifical base in phase 1.
"""
function ArtificalBase_table(b::Array{Int64,1})::SparseMatrixCSC{Int64,Int64}
    n = size(b,1)
    artificalBase_table = spzeros(n*n,1)
    for i in 0:(n-1)
        artificalBase_table[i*n + 1] = sign(b[i+1])*n
    end
    for i in 1:(n-1)
        artificalBase_table[n*n - (n - i)] = -sign(b[i])*i
    end
    return artificalBase_table
end

"""
Function: Base_table
Argument: base_list: 2 Tuple list of edges in the base.
    n: number of nodes in the network
Return: base_table: an inital look up table for the base in phase 2.
"""
function Base_table(base_list, n::Int64)::SparseMatrixCSC{Int64,Int64}
	base_table = spzeros(n*n,1)
	for (i,e) in enumerate(base_list)
            if e[1] != e[2]
            pos = (e[1] - 1)*n + 1
            while base_table[pos] != 0
                pos += 1
            end
            base_table[pos] = e[2]
            pos = (e[2] - 1)*n + 1
            while base_table[pos] != 0
                pos += 1
            end
            base_table[pos] = -e[1]
        end
    end
    return base_table
end

"""
Function: ArtificalBase_list
Argument: base_list: 2 Tuple list of edges in the base.
    n: number of nodes in the network + the artifical node
Return: artificalBase_list: an inital list of 2 Tuples for the artifical base
    in phase 1.
"""
function ArtificalBase_list(b::Array{Int64,1})::Matrix{Tuple{Int64, Int64}}
    n = size(b,1)
    artificalBase_list = Array{Tuple{Int64,Int64},2}(undef, n, 1)
    for i in 1:(n-1)
        if(sign(b[i]) == 1)
            artificalBase_list[i] = (i,n)
        else
            artificalBase_list[i] = (n,i)
        end
    end
    artificalBase_list[n] = (n,n)
    return artificalBase_list
end

"""
Function: Base_list
Argument: artificalBase_list: an inital look up table for the
    artifical base in phase 1.
    n: number of nodes in network.
Return: base_list: an inital list of 2 Tuples for the base in phase 2.
"""
function Base_list(artificalBase_list::Matrix{Tuple{Int64, Int64}}, n::Int64)::Matrix{Tuple{Int64, Int64}}
    base_list = Array{Tuple{Int64,Int64},2}(undef, n, 1)
    for (i,e) in enumerate(artificalBase_list)
        if(e != (n+1,n+1))
            if(n+1 in e)
                (e[1] != n+1) ? base_list[i] = (e[1],e[1]) : base_list[i] = (e[2],e[2])
            else
                base_list[i] = e
            end
        end
    end
    return base_list
end

"""
Function: Set_cost
Argument: base_list: an inital list of 2 Tuples for the base in phase 2.
    cost_dict: dictionary of edge costs
Return: cost: an ordered vector of cost for edges for the inital base in phase 2.
"""
function Set_cost(base_list::Matrix{Tuple{Int64, Int64}}, cost_dict::Dict{Tuple{Int64, Int64}, Int64})::Matrix{Float64}
        cost = zeros(size(base_list,1),1)
        for (i,e) in enumerate(base_list)
            if e[1]!=e[2]
                cost[i] = cost_dict[e]
            end
        end
    return cost
end

"""
Function: NonBasic_list
Argument: base_list: 2 Tuple list of edges in the base.
    edge_list: 2 Tuple list of edges in network
Return: nonBasic_list:  Tuple list of edges not in the base.
"""
function NonBasic_list(base_list::Matrix{Tuple{Int64, Int64}}, edge_list::Matrix{Tuple{Int64, Int64}})::Matrix{Tuple{Int64, Int64}}
    n = (size(edge_list,2) - size(base_list,1)) + 1             # figure out size of nb_list the + 1 is for the root edge not being counted
    nonBasic_list = Array{Tuple{Int64,Int64},2}(undef, n, 1)
    pos = 1
        for e in edge_list
            if (((e in base_list) == false)&(e[1]!=e[2]))
                println(e)
                nonBasic_list[pos] = e
                pos += 1
        end
    end
    return nonBasic_list
end

"""
Function: findPath
Argument: i, j: two nodes to find a path to
    n: number of nodes in network
    artificalBase_table: an inital look up table for the
        artifical base in phase 1.
Return: path: a 2 Tuple list of edges in path with the
    first entry being the (i,j)
"""
function Find_path(i::Int64,j::Int64,n::Int64,artificalBase_table::SparseMatrixCSC{Int64,Int64})::Array{Tuple{Int64,Int64}}
	stack = Int[]
	path = Tuple{Int64,Int64}[]
    leadingNode = Int[]
	visited = Bool[false for i in 1:n]
    prev = j
	push!(stack,-i)
	found = false
	while found == false
		#println("The Stack: $(stack)")
        #println("Visited Nodes: $(visited)")
		node = pop!(stack)
		abs_node = abs(node)
		#println("Current Node: $(abs_node)")
		if visited[abs_node] == false
			visited[abs_node] = true
			if sign(node) == -1
				#println("Adding leading node 1: ($(abs_node), $(prev))")
				push!(path,(abs_node,prev))
                push!(leadingNode,1)
            else
				#println("Adding leading node 2: ($(prev), $(abs_node))")
                push!(path,(prev,abs_node))
                push!(leadingNode,2)
            end
			prev = abs_node
			pushed = 0
			pos = (abs_node - 1)*n + 1
			neighbor = B[pos]
			while (neighbor != 0)&&(found==false)
				if abs(neighbor) == j
					found = true
					pushed += 1
                    if sign(neighbor) == 1
						#println("Found Last Arc : ($(abs_node) , $(abs(neighbor)))")
                        push!( path, (abs_node,neighbor) )
                    else
						#println("Found Last Arc: ($(neighbor) , $(abs_node))")
                        push!( path, (abs(neighbor),abs_node) )
                    end
				else
					if visited[abs(neighbor)] == false
						push!(stack,neighbor)
						pushed += 1
					end
					pos += 1
					neighbor = artificalBase_table[pos]  # Incrementing what neighbor we are looking at
				end
			end
			if(pushed == 0)  # check for dead end
				#println("$(abs_node): Dead End")
				#println("Bad Path: $(path)")
				pop!(path)
				pop!(leadingNode)
				reset = pop!(path)
				#println("Restock: $(reset)")
                #println("Leading Node: $leadingNode")
                popedLeadingNode = pop!(leadingNode)
				#println("Putting Back in Stack: $(reset[popedLeadingNode]) <= $(popedLeadingNode)")
				visited[reset[popedLeadingNode]] = false
                (popedLeadingNode == 1) ? push!(stack,-reset[popedLeadingNode]) : push!(stack,reset[popedLeadingNode])
				#println((popedLeadingNode == 1) ? 2 : 1)
				prev = reset[(popedLeadingNode == 1) ? 2 : 1]
			end
		end
	end
	return path
end

"""
Function: BaseMatrix
Argument: edge_list: 2 Tuple list of edges for either phase 1 or phase 2.
    n: number of nodes in the network
Return: B_matrix: sparse adjacency matrix
"""
function BaseMatrix(edge_list::Matrix{Tuple{Int64, Int64}},n::Int64)::SparseMatrixCSC{Int64,Int64}
	B_matrix = spzeros(n,n)
	for (k,e) in enumerate(edge_list)
		B_matrix[e[1],k] = 1
		B_matrix[e[2],k] = -1
	end
	return B_matrix
end

"""
Function: Update_BaseMatrix!
Argument: B_matrix: sparse adjacency matrix. edge_list: 2 Tuple list of edges for either phase 1 or phase 2.
    exitingEdge, enteringEdge: 2 Tuple of the edge exiting/entering
Return: None
"""
function Update_BaseMatrix!(B_matrix::SparseMatrixCSC{Int64,Int64},edge_list::Matrix{Tuple{Int64, Int64}},exitingEdge::Tuple{Int64, Int64},enteringEdge::Tuple{Int64, Int64})
    i = 1
    while B_edgeList[i] != exitingEdge
        i += 1
    end
    B_edgeList[i] = enteringEdge
    B_matrix[exitingEdge[1],i] = 0
    B_matrix[exitingEdge[2],i] = 0

    B_matrix[enteringEdge[1],i] = 1
    B_matrix[enteringEdge[2],i] = -1
end

"""
Function: Get_direction
Argument: path: a 2 Tuple list of edges in path with the
    first entry being the (i,j)
Return: forward: this bool is to determin if the path gose around
 the cycle forwrds or backwards
"""
function Get_direction(path::Array{Tuple{Int64,Int64}})::Bool
	if (path[1][2]==path[2][2])||(path[1][2]==path[2][1])
		forward = true
	else
		forward = false
	end
	return forward
end

"""
Function: Get_direction
Argument: path: a 2 Tuple list of edges in path with the
    first entry being the (i,j)
    n: number of no
Return: flow_list: a list with +1 if the edge is in the same direction
    of the new edge -1 otherwise.
    valid: if flow list has all +1 then the algorithm can not proceed
"""
function Get_flow(path::Array{Tuple{Int64,Int64}})
	valid = false;
	prev = path[1]
	withFlow = true
	flow_list = zeros(Int, size(path,1))
	flow_list[1] = 1
	#determin the ordering of the path
	forward = Get_direction(path)
	#println("What Direction Are we Going??? $(forward)")
	if forward
		range = (2:size(path, 1))
	else
		range = (size(path, 1):(-1):2)
	end
	i=2; j=1
	for k in range
		if ( (prev[i] != path[k][j]) && (withFlow == true))
			#println("$(prev) --> $(path[k]) if")
			flow_list[k] = -1
			withFlow = false
			valid = true
			i=1; j=2
		elseif (prev[i] == path[k][j])&&(withFlow == false)
			#println("$(prev) --> $(path[k]) elif")
			flow_list[k] = -1
		else
			#println("$(prev) --> $(path[k]) else")
            flow_list[k] = 1
			withFlow = true
			i=2; j=1
		end
		prev = path[k]
	end
    return flow_list, valid
end

"""
Function: ReducedCost
Argument: w: dual values. nonBasic_list:  Tuple list of edges not in the base.
    phase: what phase is this function being called from 1 or 2.
    cost: for phase 1 send 0 for phase 2 pass the cost_dict
Return: max: largest reduced cost
    newBasicEdge: edge to be added to base
"""
function ReducedCost(w, nonBasic_list, phase, cost)
    max = typemin(Float64)
    newBasicEdge = (0,0)
    c = 0.0
    for (i,e) in enumerate(nonBasic_list)
        if e != (0,0)                               #(0,0) for deleted nb_edges
            (phase == 1) ? c = 0 : c = cost[e]
            rij = w[e[1]] - w[e[2]] - c
			#println("$(e): $(rij) =  $(w[e[1]] - w[e[2]] - c)")
            if(rij > max)
                max = rij
                newBasicEdge = e
            end
        end
    end
    return max, newBasicEdge
end

"""
Function: Get_Indexes
Argument:
Return:
"""
function Get_Indexes(path,B_edgeList)
    indexes = Int[]
    for i in 2:size(path,1)
        for j in 1:size(B_edgeList,1)
            if(path[i] == B_edgeList[j])
                push!(indexes,j)
                break
            end
        end
    end
    return indexes
end

"""
Function:
Argument:
Return:
"""
function getDeleta(flow, path, indexes, x)
    #Find the min delta i.e smallest shipping pattern
    delta = typemax(Int64)
    min_edge = (0,0)
    min_index = 0

	#determin the ordering of the path
	forward = Get_direction(path)
	if forward
		range = (2:size(path, 1))
	else
		range = (size(path, 1):(-1):2)
	end

    for i in range
		#println("Flow on $(path[i]): $(x[i-1])")
        if (flow[i] == -1)
            k = indexes[i-1]
            if x[k] < delta
                delta = x[k]
                min_edge = path[i]
                min_index = k
            end
        end
    end
    return delta, min_edge, min_index
end

"""
Function:
Argument:
Return:
"""
function updateFlowPattern!(flow, path, indexes, min_index, delta, x)
    # First update the incomming edged spot in outgoing spot
    x[min_index] = delta
    for i in 2:size(flow,1)
		k = indexes[i-1]
        if k != min_index
            #println("Index: $(i) and Flow: $(flow[i])")
			if flow[i] == 1
				temp_x = x[k]
				x[k] += delta
				#println("$(temp_x) ---> $(x[i])")
			else
				temp_x = x[k]
				x[k] -= delta
				#println("$(temp_x) ---> $(x[i])")
			end
        end
    end
end

"""
Function:
Argument:
Return:
"""
function RemoveFromTable!(i,j,n,B)
    tempStack = Int64[]
    pos = (i - 1)*n + 1
	#println("pos: $(pos)")
	for k in 1:n
		if abs(B[pos]) ==j
			#println("Found abs(B[pos]) at pos: $(pos)")
			break
		end
		pos += 1
    end
    B[pos] = 0
	#println(B)
    keep_pos = pos + 1
    while(B[keep_pos] != 0)
		#println("Pushing $(B[keep_pos]) from pos: $(keep_pos)")
        push!(tempStack, B[keep_pos])
		B[keep_pos] = 0
        keep_pos +=1
    end

    if(size(tempStack,1) > 0)
        for k in 1:size(tempStack,1)
			#println("Placing $(tempStack[k]) back on table at pos: $(pos)")
            B[pos] = tempStack[k]
            pos+=1
        end
	end
end

"""
Function:
Argument:
Return:
"""
function AddToTable!(u, v,node_sign,B)
	pos = (u - 1)*n + 1
	while(B[pos] != 0)
		pos += 1
	end
	B[pos] = node_sign*v
end
