module CustomLDA
  export LDA

  const CONCENTRATION = 0.5

  function LDA(probabilities_array, counts=ones(Int32, length(probabilities_array)))
    tag_count = length(probabilities_array)

    prior = Array{Float64}(tag_count) # aka theta
    posterior= Array{Float64}(tag_count)
    for i in 1:tag_count
      prior[i] = 1.0/tag_count
      posterior[i] = 0.0
    end
    iterations = 100
    temp_norm_total = 0.0
    for iteration in 1:iterations
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
        temp_norm_total += posterior[real_index] + CONCENTRATION
      end
      for real_index in 1:tag_count
        prior[real_index] = (posterior[real_index] + CONCENTRATION) / temp_norm_total
        posterior[real_index] = 0.0
      end
    end

    # Return value
    most_likely_real_for_each_obs = Array{Tuple{Int32, Float64}}(tag_count)
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
