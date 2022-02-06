type t [@@deriving sexp]
val to_float        : ?idx:int -> t -> float
val to_float_array  : t -> float array
val of_float        : float -> t
val of_float_array  : dim:int -> float array -> t
val to_string       : t -> string
val to_bytes        : t -> bytes
val of_bytes        : bytes -> t
val dimension       : t -> int

