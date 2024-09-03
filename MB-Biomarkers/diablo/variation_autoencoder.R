library(keras)

# Define the SMILES strings for your drug molecules
group3_drugs <- fread(file.path(output_dir, "group3_finalcids_prop.csv"))
smiles_data <- group3_drugs$isosmiles  

# Model Building
latent_dim <- 128
max_len <- 100  
charset_length <- 20 

# Encoder
encoder_input <- layer_input(shape = c(max_len, charset_length))
encoder_lstm <- layer_lstm(units = 256)(encoder_input)
z_mean <- layer_dense(units = latent_dim)(encoder_lstm)
z_log_var <- layer_dense(units = latent_dim)(encoder_lstm)

# Sampling layer
sampling <- function(args) {
  z_mean <- args[[1]]
  z_log_var <- args[[2]]
  epsilon <- backend$random_normal(shape = k_shape(z_mean),
                                   mean = 0, stddev = 1)
  z <- backend$mean(z_mean + backend$exp(0.5 * z_log_var) * epsilon)
  return(z)
}

z <- layer_lambda(function(x) sampling(x))(list(z_mean, z_log_var))


# Decoder
decoder_input <- layer_input(shape = c(latent_dim))
decoder_lstm <- layer_lstm(units = 256, return_sequences = TRUE)(decoder_input)
decoder_output <- layer_time_distributed(layer_dense(units = charset_length,
                                                     activation = "softmax"))(decoder_lstm)

# Compile the VAE model
vae <- keras_model(inputs = encoder_input, outputs = decoder_output)
kl_loss <- -0.5 * k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var))
vae_loss <- k_mean(k_binary_crossentropy(encoder_input, decoder_output)) + kl_loss
compile(vae, loss = "binary_crossentropy", optimizer = "adam")

# Model Training
history <- fit(vae, x = smiles_data, y = smiles_data,
               epochs = 100, batch_size = 128)

# Generating New Molecules
# Generating New Molecules
# Function to generate molecules from latent space
generate_molecules <- function(vae_model, latent_dim, charset_length, max_len, n_samples = 10) {
  # Generate random samples from the latent space
  latent_samples <- array(runif(n_samples * latent_dim, min = -1, max = 1), dim = c(n_samples, latent_dim))
  
  # Decode the samples to generate molecules
  decoder_input_layer <- layer_input(shape = c(latent_dim))
  decoder_lstm_layer <- layer_lstm(units = 256, return_sequences = TRUE)
  decoder_output_layer <- layer_time_distributed(layer_dense(units = charset_length,
                                                             activation = "softmax"))
  decoder_output <- decoder_output_layer(decoder_lstm_layer(decoder_input_layer))
  decoder_model <- keras_model(inputs = decoder_input_layer, outputs = decoder_output)
  decoder_model$set_weights(vae_model$get_layer(index = 4)$get_weights()) # Assuming decoder is the fourth layer
  
  generated_molecules <- vector(mode = "list", length = n_samples)
  for (i in 1:n_samples) {
    latent_sample <- matrix(latent_samples[i, ], nrow = 1)
    decoded_molecule <- predict(decoder_model, latent_sample)
    generated_molecules[[i]] <- decode_molecule(decoded_molecule, charset_length)
  }
  
  return(generated_molecules)
}

# Function to decode one-hot encoded molecule representation to SMILES string
decode_molecule <- function(molecule, charset_length) {
  decoded_molecule <- ""
  for (j in 1:dim(molecule)[2]) {
    char_probs <- molecule[1, j, ]
    char_idx <- sample.int(charset_length, size = 1, prob = char_probs)
    decoded_molecule <- paste(decoded_molecule, charset[char_idx], sep = "")
  }
  return(decoded_molecule)
}

# Assuming charset and max_len are defined
generated_molecules <- generate_molecules(vae, latent_dim, charset_length, max_len, n_samples = 10)
print(generated_molecules)

