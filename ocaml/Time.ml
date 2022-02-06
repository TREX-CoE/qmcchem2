let of_sec s =
  Unix.gmtime s

let to_sec t =
  let sec  = t.Unix.tm_sec
  and min  = t.Unix.tm_min
  and hour = t.Unix.tm_hour
  and mday = t.Unix.tm_mday
  in
  sec +
  min  * 60 +
  hour * 60 * 60 + 
  (mday-1) * 60 * 60 * 24 

let string_of_t t = 
  let mday = t.Unix.tm_mday - 1 in
  let sec  = t.Unix.tm_sec
  and min  = t.Unix.tm_min
  and hour = t.Unix.tm_hour + 24*mday
  in
  Printf.sprintf "%2d:%2.2d:%2.2d" hour min sec

let string_of_date t = 
  let year = 1900 + t.Unix.tm_year in
  let mon  = t.Unix.tm_mon  in
  let mday = t.Unix.tm_mday in
  let sec  = t.Unix.tm_sec
  and min  = t.Unix.tm_min
  and hour = t.Unix.tm_hour 
  in
  let month = 
    match mon with 
      | 0  -> "Jan" | 1  -> "Feb" | 2  -> "Mar" | 3  -> "Apr"
      | 4  -> "May" | 5  -> "Jun" | 6  -> "Jul" | 7  -> "Aug"
      | 8  -> "Sep" | 9  -> "Oct" | 10 -> "Nov" | 11 -> "Dec"
      | _ -> assert false
  in
  Printf.sprintf "%2d %3s %4d - %2d:%2.2d:%2.2d" mday month year hour min sec


let string_of_now () = 
  Unix.gettimeofday ()
  |> Unix.localtime 
  |> string_of_date


let string_of_sec s = 
  of_sec s 
  |> string_of_t 

