pub fn get_rng() -> impl ark_std::rand::Rng {
    #[cfg(test)]
    return ark_std::test_rng();

    #[cfg(not(test))]
    return ark_std::rand::thread_rng();
}
