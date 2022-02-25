import tensorflow as tf 
from tensorflow.keras.layers import Layer  

class CombineCfacLayer(Layer):
   """
   Creates the extra layer to fit SMEFT c-factors
   on top of NNPDF4.0's architecture.  
   
   Parameters
   ----------
            ncfacs: int 
               Defines the number of Wilson 
               coefficients to be included in a the fit. 
   """
   
   def __init__(self, ncfacs):
      super().__init__()

      init_value = tf.random_normal_initializer()
      self.w = tf.Variable(
         initial_value=init_value(shape=(ncfacs,), dtype='float32'),
         trainable = True
         )
   
   def call(self, inputs, cfactor_values):
      return (1 + tf.reduce_sum(w[:, tf.newaxis] * cfactor_values, axis=0)) * inputs

                              




