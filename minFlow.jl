

"""
Function: ArtificalBase_table
Argument: b a vector of node supply and demand
Return: B a lookup hash table
"""
function ArtificalBase_table(b::Array{Int64,1})::SparseMatrixCSC{Int64,Int64}
    n = size(b,1)
    B = spzeros(n*n,1)
    for i in 0:(n-1)
        B[i*n + 1] = sign(b[i+1])*n
    end
    for i in 1:(n-1)
        B[n*n - (n - i)] = -sign(b[i])*i
    end
    return B
end

"""


"""
function Base_table(base_list, n)::SparseMatrixCSC{Int64,Int64}
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
Function:
"""
function ArtificalBase_list(b::Array{Int64,1})::Array{Tuple{Int64,Int64},1}
    n = size(b,1)
    artificalBase_list = Tuple{Int64,Int64}[]
    for i in 1:(n-1)
        if(sign(b[i]) == 1)
        push!(artificalBase_list, (i,n))
        else
            push!(artificalBase_list, (n,i))
        end
    end
    push!(artificalBase_list, (n,n))
    return artificalBase_list
end



"""
Function:
"""
function Base_list(artificalBase_list::Array{Tuple{Int64,Int64},1})::Array{Tuple{Int64,Int64},1}
    base_list = Tuple{Int64,Int64}[]
    for (i,e) in enumerate(artificalBase_list)
        if(e != (n,n))
            if(n in e)
                (e[1] != n) ? push!(base_list, (e[1],e[1])) : push!(base_list, (e[2],e[2]))
            else
                push!(base_list, e)
            end
        end
    end
    return base_list
end

"""
Function:
"""
function Set_cost(base_list, cost_dict)
        cost = zeros(size(base_list,1),1)
        for (i,e) in enumerate(base_list)
            if e[1]!=e[2]
                cost[i] = cost_dict[e]
            end
        end
    return cost
end

"""
Function:
"""
function NonBasic_list(base_list, edge_list)
    nonBasic_list = Tuple{Int64,Int64}[]
        for e in edge_list
            if ((e in base_list) == false)
                push!(nonBasic_list, e)
        end
    end
    return nonBasic_list
end


"""
Function: findPath
Find a path from node i to node j if B is a ArtificalBase_table
"""
function findPath(i::Int64,j::Int64,n::Int64,B::SparseMatrixCSC{Int64,Int64})::Array{Tuple{Int64,Int64}}

	stack = Int[]
	path = Tuple{Int64,Int64}[]
    leadingNode = Int[]
	visited = Bool[false for i in 1:n]

    prev = j
	push!(stack,-i)
	found = false
	while found == false
		#println("The Stack: $(stack)")
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
					neighbor = B[pos]  # Incrementing what neighbor we are looking at
				end
			end

			if(pushed == 0)  # check for dead end
				#println("$(abs_node): Dead End")
				#println("Bad Path: $(path)")
				pop!(path)
				pop!(leadingNode)
				reset = pop!(path)
				#println("Restock: $(reset)")
                popedLeadingNode = pop!(leadingNode)
				#println("Putting Back in Stack: $(reset[popedLeadingNode]) <= $(popedLeadingNode)")
				visited[reset[popedLeadingNode]] = false
				push!(stack,reset[popedLeadingNode])
				#println((popedLeadingNode == 1) ? 2 : 1)
				prev = reset[(popedLeadingNode == 1) ? 2 : 1]
			end

		end
	end

	return path

end

"""


"""
function Base_Matrix(edgeList::Array{Tuple{Int64,Int64},1},n::Int64)::SparseMatrixCSC{Int64,Int64}
	B_matrix = spzeros(n,n)
	for (k,e) in enumerate(edgeList)
		B_matrix[e[1],k] = 1
		B_matrix[e[2],k] = -1
	end
	return B_matrix
end

"""


"""
function getDirection(path)
	if (path[1][2]==path[2][2])||(path[1][2]==path[2][1])
		forward = true
	else
		forward = false
	end
	return forward
end

"""


"""
function getFlow(path,n)
	valid = false;
	prev = path[1]
	withFlow = true
	D = zeros(Int, size(path,1))
	D[1] = 1
	#determin the ordering of the path
	forward = getDirection(path)
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
			D[k] = -1
			withFlow = false
			valid = true
			i=1; j=2
		elseif (prev[i] == path[k][j])&&(withFlow == false)
			#println("$(prev) --> $(path[k]) elif")
			D[k] = -1
		else
			#println("$(prev) --> $(path[k]) else")
            D[k] = 1
			withFlow = true
			i=2; j=1
		end
		prev = path[k]
	end
    return D, valid
end

"""


"""
function reducedCost(w, NB_edgeList, Cost)
    max = typemin(Float64)
    newBasicEdge = (0,0)
    for (i,e) in enumerate(NB_edgeList)
        if e != (0,0)
            rij = w[e[1]] - w[e[2]] - C[e[1],e[2]]
			#println("$(e): $(rij) =  $(w[e[1]]) - $(w[e[2]]) - $(C[e[1],e[2]])")
            if(rij > max)
                max = rij
                newBasicEdge = e
            end
        end
    end
    return max, newBasicEdge
end

"""


"""
function getIndexes(path,B_edgeList)
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


"""
function getDeleta(flow, path, indexes, x)
    #Find the min delta i.e smallest shipping pattern
    delta = typemax(Int64)
    min_edge = (0,0)
    min_index = 0

	#determin the ordering of the path
	forward = getDirection(path)
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


"""
function AddToTable!(u, v,node_sign,B)
	pos = (u - 1)*n + 1
	while(B[pos] != 0)
		pos += 1
	end
	B[pos] = node_sign*v
end


"""
# find the exiting edge in the list
# update and use the index to update the Matrix
"""
function updateBaseMatrix!(B_matrix,B_edgeList,exitingEdge,enteringEdge)
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
