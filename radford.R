radford = function (UgradU, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q), 0, 1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  ugc = UgradU(q)
  p = p - epsilon * ugc$dudq / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    ug = UgradU(q)
    
    cK = sum(p^2) / 2
    cat(ug$u + cK, "|", ug$u, ",", cK, "|", q, "\n")
    
    #print(paste(q, p, ug$u, ugc$u, sum(p^2) / 2, sum(current_p^2) / 2))
    if (i!=L) {
      p = p - epsilon * ug$dudq
    }
  }
  # Make a half step for momentum at the end.
  ug = UgradU(q)
  p = p - epsilon * ug$dudq / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = ugc$u
  current_K = sum(current_p^2) / 2
  proposed_U = ug$u
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  #print(current_U + current_K)
  #print(proposed_U + proposed_K)
  #print(proposed_K)
  #if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  #{
  return (q) # accept
  #}
  #else
  #{
  #  return (current_q) # reject
  #}
}