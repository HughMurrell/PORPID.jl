module CustomLDA
  export LDA

  const DEFAULT_CONCENTRATION = 0.5
  const EPSILON = 0.00000000000000001

  function LDA(probabilities_array, counts=ones(Int32, length(probabilities_array)); concentration=DEFAULT_CONCENTRATION)
    tag_count = length(probabilities_array)
    converged_epsilon =  EPSILON / tag_count

    prior = Vector{Float64}(undef, tag_count) # aka theta
    posterior = zeros(Float64, tag_count)
    for i in 1:tag_count
      prior[i] = 1.0/tag_count
    end
    max_iterations = 1000
    temp_norm_total = 0.0
    for iteration in 1:max_iterations
      # φ = Map[norm[θ*#] &, CL]
      for obs in 1:tag_count
        # Calculating the total for this row
        temp_norm_total = 0.0
        for (real_index, prob) in probabilities_array[obs]
          temp_norm_total += prior[real_index] * prob
        end
        normalized_count = counts[obs] / temp_norm_total
        for (real_index, prob) in probabilities_array[obs]
          posterior[real_index] += (prior[real_index] * prob) * normalized_count
        end
      end
      temp_norm_total = 0.0
      for real_index in 1:tag_count
        temp_norm_total += posterior[real_index] + concentration
      end
      converged = true
      for real_index in 1:tag_count
        new_prior = (posterior[real_index] + concentration) / temp_norm_total
        converged = (converged && abs(new_prior-prior[real_index]) < converged_epsilon)
        prior[real_index] = new_prior
        posterior[real_index] = 0.0
      end
      if (converged)
        println(stderr, "Converged after $(iteration) iterations")
        break
      end
    end

    # Return value
    most_likely_real_for_each_obs = Array{Tuple{Int32, Float64}}(undef, tag_count)
    for obs in 1:tag_count
      best_prob = 0.0
      best_index = -1
      temp_norm_total = 0.0
      for (real_index, prob) in probabilities_array[obs]
        posterior = prior[real_index] * prob
        temp_norm_total += posterior
        if posterior > best_prob
          best_prob = posterior
          best_index = real_index
        end
      end
      most_likely_real_for_each_obs[obs] = (best_index, best_prob / temp_norm_total)
    end
    return most_likely_real_for_each_obs
  end

  function norm(v)
    return v ./ sum(v)
  end

end
