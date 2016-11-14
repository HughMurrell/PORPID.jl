module SparseICMapLDA
  export LDA

  const CONCENTRATION = 0.5

  function LDA(cl_matrix, counts=zeros(Int32, 0, 1))
    rows = size(cl_matrix, 1)
    cols = size(cl_matrix, 2)
    #theta must be a 1-by-x matrix to make element-wise multiplication work correctly
    theta = ones(Float32, 1, cols) / cols
    phi = spzeros(Float32, rows, cols)
    #counts must be a <rows>-by-1 matrix for the same reason
    if size(counts, 1) == 0
      counts = ones(Int32, rows, 1)
    end

    iterations = 100
    for i in 1:iterations
      for r in 1:rows
        phi[r,:] = norm(theta .* cl_matrix[r,:])
      end
      for c in 1:cols
        theta[1,c] = sum(phi[:,c] .* counts) + CONCENTRATION
      end
      theta = norm(theta)
      println(STDERR, i)
    end

    return phi
  end

  function norm(v)
    return v ./ sum(v)
  end

end
