import hazel

iterator = hazel.Iterator(use_mpi=True)
mod = hazel.Model('conf_synth_parallel.ini', working_mode='synthesis')
iterator.use_model(model=mod)
iterator.run_all_pixels()
