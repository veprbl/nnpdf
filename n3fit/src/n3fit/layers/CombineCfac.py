import tensorflow as tf 
from tensorflow.keras.layers import Layer  

class CombineCfacLayer(Layer):
   """
   Creates SIMUnet's combination layer.  
   """
   
   def __init__(self, ncfacs, cfac_units, cfac_labels, quad_labels):
      # Initialize a Layer instance
      super().__init__()    
      # Initialize vector of trainable weights
      self.w = tf.Variable(
         initial_value=tf.zeros(shape=(ncfacs,), dtype='float32'),
         trainable = True
         )

      self.cfac_labels = cfac_labels
      self.quad_labels = quad_labels
      self.quad_indices = [self.cfac_labels.index(i) for i in self.quad_labels]
      self.scale = cfac_units
   
   def __call__(self, inputs, cfactor_values, quad_cfactor_values):
      """
      The operation that the SIMUnet layer performs on inputs.
      At this point we only include linear interference with the SMEFT. 
      Parameters
      ----------
         inputs: float 
            Represents the SM (c_factor_values = 0)
            theoretical prediction for a given observable. 
         cfactor_values: tf.Variable
            Set of trainable SMEFT C-factors that affect a given
            observable. 

      """
      linear = tf.reduce_sum(self.w[:, tf.newaxis] * cfactor_values, axis=0) / self.scale
      if self.quad_indices:
         # https://stackoverflow.com/questions/46881006/slicing-tensor-with-list-tensorflow/51139591
         square_weights = tf.stack([self.w[i] for i in self.quad_indices])**2
         quadratic = tf.reduce_sum(square_weights[:, tf.newaxis] * quad_cfactor_values, axis=0) / (self.scale)**2
      else:
         quadratic = tf.zeros(linear)
      return (1 + linear + quadratic) * inputs

                              




