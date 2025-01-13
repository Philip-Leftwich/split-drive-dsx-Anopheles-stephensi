source("scripts/fecundity_analysis.R")
library(patchwork)

data %>% 
  group_by(Cross) %>% 
  summarise(sum = sum(`Larvae number`, na.rm = T)/50) %>% 
  mutate(relative_fitness = c(35.1/63.8, 12/63.8, 10.7/53.5, 53.5/53.5, 63.8/63.8)) %>% 
  mutate(fitness_cost = 1-relative_fitness)

# Homing Endonuclease Gene (HEG) Spread Model with Sex-Specific Parameters====

## 1. Define Parameters====

# Fitness costs for homozygotes
s_m <- 0.812   # Fitness cost for male homozygotes (HEG/HEG)
s_f <- 1   # Fitness cost for female homozygotes (HEG/HEG)

# Fitness costs for heterozygotes
sh_m <- 0.45  # Fitness cost for male heterozygotes (HEG/Wild)
sh_f <- 0.8 # Fitness cost for female heterozygotes (HEG/Wild)

# Homing rates
e_m <- 0.98   # Homing rate in males
e_f <- 0.98   # Homing rate in females


# 2. Initialize Frequencies====

q <- .1

# Initial HEG frequencies
q_m <- q   # Initial frequency in males
q_f <- q   # Initial frequency in females


# 3. Set Convergence Criteria====


threshold <- 1e-6       # Convergence threshold
max_generations <- 100 # Maximum number of generations to prevent infinite loops
generation <- 0          # Generation counter


# 4. Iterative Process====


# Create vectors to store frequencies over generations (optional for analysis)
freq_m <- c(q_m)
freq_f <- c(q_f)
loads <- c()


# Iterate until convergence or until maximum generations are reached
while (generation < max_generations) {
  
  # Calculate allele frequencies
  p_m <- 1 - q_m  # Frequency of wild-type allele in males
  p_f <- 1 - q_f  # Frequency of wild-type allele in females
  
  # ---------------------------
  # Update for Males
  # ---------------------------
  
  # Numerator for males
  numerator_m <- (1 - s_m) * q_m*q_f + (1 - sh_m)* 0.5 * ((p_m*q_f)+(p_f*q_f)) * (1 + e_m)
  
  # Denominator for males
  denominator_m <- 1 - s_m * q_m*q_f - sh_m * ((p_m*q_f)+(p_f*q_f))
  
  # New HEG frequency in males
  ql_m <- numerator_m / denominator_m
  
  # ---------------------------
  # Update for Females
  # ---------------------------
  
  # Numerator for females
  numerator_f <- (1 - s_f) * q_m*q_f + (1 - sh_f)* 0.5 * ((p_m*q_f)+(p_f*q_f)) * (1 + e_f)
  
  # Denominator for females
  denominator_f <- 1 - s_f * q_m*q_f - sh_f * ((p_m*q_f)+(p_f*q_f))
  
  # New HEG frequency in females
  ql_f <- numerator_f / denominator_f
  
  # ---------------------------
  # Check for Convergence
  # ---------------------------
  
  # Calculate changes
  delta_m <- abs(ql_m - q_m)
  delta_f <- abs(ql_f - q_f)
  
  q <- (ql_m+ql_f)/2
  
  # Update frequencies
  q_m <- ql_m
  q_f <- ql_f
  
  # Store frequencies and load (optional)
  freq_m <- c(freq_m, q_m)
  freq_f <- c(freq_f, q_f)
  load_current <- (s_m * q_m) + (s_f * q_f)
  loads <- c(loads, load_current)
  
  # Increment generation counter
  generation <- generation + 1
  
  # Check if both changes are below the threshold
  if (delta_m < threshold && delta_f < threshold) {
    break
  }
}


# 5. Output Results=====


cat("Converged after", generation, "generations:\n")
cat(sprintf("HEG frequency in males: %.4f\n", q_m))
cat(sprintf("HEG frequency in females: %.4f\n", q_f))

# Calculate and display total population load
total_load <- (s_m * q_m) + (s_f * q_f)
cat(sprintf("Total population load: %.4f (%.2f%%)\n", total_load, total_load * 100))


# 6. (Optional) Plotting Results====


data <- tibble(freq_m, freq_f, generation = (0:generation))

plot_1 <- data %>% 
  pivot_longer(cols = freq_f:freq_m, names_to = "sex", values_to = "frequency") %>% 
  ggplot(aes(x = generation, y = frequency, linetype = sex))+
  geom_line()+
  scale_y_continuous(limits = c(0,1), n.breaks = 10)+
  theme_bw()



# Finding acceptable het fitness loads====

sh <- seq(0,1, by =.05)


data_list <-  vector("list", length = length(sh))

for (i in 1:length(sh)){
  q <- .1
  
  # Initial HEG frequencies
  q_m <- q   # Initial frequency in males
  q_f <- q   # Initial frequency in females
  
  # Create vectors to store frequencies over generations (optional for analysis)
  freq_m <- c(q_m)
  freq_f <- c(q_f)
  generation <- 0
  generations <- c(generation)
  loads <- c()
  sh_1 <- sh[[i]]
  sh_m <- sh[[i]]
  sh_f <- sh[[i]]
  
  # Iterate until convergence or until maximum generations are reached
  while (generation < max_generations) {
    
    # Calculate allele frequencies
    p_m <- 1 - q_m  # Frequency of wild-type allele in males
    p_f <- 1 - q_f  # Frequency of wild-type allele in females
    
    # ---------------------------
    # Update for Males
    # ---------------------------
    
    # Numerator for males
    numerator_m <- (1 - s_m) * q_m*q_f + (1 - sh_m)* 0.5 * ((p_m*q_f)+(p_f*q_f)) * (1 + e_m)
    
    # Denominator for males
    denominator_m <- 1 - s_m * q_m*q_f - sh_m * ((p_m*q_f)+(p_f*q_f))
    
    # New HEG frequency in males
    ql_m <- numerator_m / denominator_m
    
    # ---------------------------
    # Update for Females
    # ---------------------------
    
    # Numerator for females
    numerator_f <- (1 - s_f) * q_m*q_f + (1 - sh_f)* 0.5 * ((p_m*q_f)+(p_f*q_f)) * (1 + e_f)
    
    # Denominator for females
    denominator_f <- 1 - s_f * q_m*q_f - sh_f * ((p_m*q_f)+(p_f*q_f))
    
    # New HEG frequency in females
    ql_f <- numerator_f / denominator_f
    
    # ---------------------------
    # Check for Convergence
    # ---------------------------
    
    # Calculate changes
    delta_m <- abs(ql_m - q_m)
    delta_f <- abs(ql_f - q_f)
    
    q <- (ql_m+ql_f)/2
    
    # Update frequencies
    q_m <- ql_m
    q_f <- ql_f
    
    # Store frequencies and load (optional)
    freq_m <- c(freq_m, q_m)
    freq_f <- c(freq_f, q_f)
    load_current <- (s_m * q_m) + (s_f * q_f)
    loads <- c(loads, load_current)
    
    # Increment generation counter
    
    generations <- c(generations, generation)
    generation <- generation + 1
    data <- tibble(freq_m, freq_f, sh_1, generations)
    data_list[[i]] <- data
  }
}


sim_data <- bind_rows(data_list) %>% pivot_longer(cols = freq_m:freq_f, names_to = "sex", values_to = "freq")


plot_2 <- sim_data %>% 
  ggplot(aes(x = generations,
             y = freq))+
  geom_line(aes(group = interaction(sh_1,sex),
                linetype = sex))+
  facet_wrap(~sh_1)+
  theme_bw()

# Combine plots====

plot_1 +plot_2 + plot_layout(guides = 'collect')
