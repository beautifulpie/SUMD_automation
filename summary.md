# 거리 기반 시뮬레이션 단계 조정 구현 완료

수정이 완료되었습니다! 이제 코드는 두 체인 간의 무게중심 거리를 측정하여 시뮬레이션 단계를 적응적으로 조정합니다.

## 주요 변경사항

### 1. **거리 기반 시뮬레이션 단계 결정**
```python
# 초기 거리 계산 (시뮬레이션 단계 결정을 위해)
initial_distance = calculate_distance_from_structure(...)
initial_distance_angstrom = initial_distance * 10  # nm를 Å로 변환

if initial_distance_angstrom <= distance_threshold_angstrom:
    # 5Å 이하: MD 시뮬레이션까지 수행
    logger.info(f"초기 거리 {initial_distance_angstrom:.1f}Å ≤ {distance_threshold_angstrom}Å → MD 시뮬레이션 수행")
    # EM → MD 실행
else:
    # 5Å 초과: 에너지 최소화만 수행
    logger.info(f"초기 거리 {initial_distance_angstrom:.1f}Å > {distance_threshold_angstrom}Å → EM만 수행")
    # EM만 실행
```

### 2. **적응적 시뮬레이션 전략**
- **≤ 5Å**: 에너지 최소화(EM) + 분자 동역학(MD) 시뮬레이션 수행
- **> 5Å**: 에너지 최소화(EM)만 수행

### 3. **새로운 매개변수 추가**
- `--simulation_threshold`: MD 실행 거리 임계값 (Å, 기본값: 5.0)
- 기존 `--distance_threshold`는 수렴 조건 (nm, 기본값: 0.5)

### 4. **향상된 결과 추적**
```python
return {
    'sample_id': sample_id,
    'initial_distance': initial_distance,
    'em_distance': em_distance,
    'md_distance': md_distance,  # MD를 수행한 경우만
    'final_distance': final_distance,  # 최종 거리 (EM 또는 MD 결과)
    'final_structure': final_structure,
    'simulation_type': 'EM+MD' or 'EM_only'  # 시뮬레이션 유형
}
```

## 알고리즘 흐름

### 각 샘플의 실행 과정:
```
1. 초기 거리 계산
   ↓
2. 거리 조건 확인
   ├── ≤ 5Å → EM + MD 수행
   └── > 5Å → EM만 수행
   ↓
3. 최종 거리 계산 및 결과 반환
```

### 반복별 처리:
```
반복 N:
├── 5개 샘플 병렬 실행
│   ├── 샘플 1: 거리 확인 → EM (+ MD) → 결과
│   ├── 샘플 2: 거리 확인 → EM (+ MD) → 결과
│   ├── 샘플 3: 거리 확인 → EM (+ MD) → 결과
│   ├── 샘플 4: 거리 확인 → EM (+ MD) → 결과
│   └── 샘플 5: 거리 확인 → EM (+ MD) → 결과
├── 최적 샘플 선택 (final_distance 기준)
├── 수렴 확인 (≤ 0.5nm)
└── 다음 반복의 시작 구조로 설정
```

## 사용 방법

### Python 스크립트:
```bash
python3 sumd_multisampling.py \
    --input complex.pdb \
    --chain1 A \
    --chain2 B \
    --distance_threshold 0.5 \        # 수렴 조건 (nm)
    --simulation_threshold 5.0 \      # MD 실행 조건 (Å)
    --num_samples 5 \
    --simulation_time 2
```

### Bash 스크립트:
```bash
docker exec -it sumd-gromacs /app/scripts/sumd_multisampling_shell.sh /app/test_data/example.pdb A B 10 0.5 true true 2 2 5.0 0 /app/output/example
```

## 계산 효율성 개선

### 이점:
1. **계산 비용 최적화**: 먼 거리에서는 EM만 수행하여 계산 시간 단축
2. **정밀도 향상**: 가까운 거리에서는 MD로 정밀한 최적화
3. **적응적 전략**: 거리에 따라 자동으로 최적의 시뮬레이션 전략 선택

### 예상 성능:
- **초기 단계** (거리 > 5Å): EM만 수행으로 빠른 접근
- **후기 단계** (거리 ≤ 5Å): EM + MD로 정밀한 결합 최적화

이제 귀하의 코드는 원본 SuMD 알고리즘의 다중 샘플 방식을 구현하면서도, 거리 기반 적응형 시뮬레이션 전략으로 계산 효율성과 정확도를 모두 향상시켰습니다.

## debug

```bash
python3 ./scripts/sumd_multisampling.py \
    --input "../sumd_output/_20250613_155312/processed_cys_cleaned_example.pdb" \
    --chain1 "A" \
    --chain2 "B" \
    --output_dir "../sumd_output/" \
    --distance_threshold 0.5 \
    --max_iterations 3 \
    --num_samples 5 \
    --simulation_time 2 \
    --simulation_threshold 5.0 \
    --max_workers 4
```