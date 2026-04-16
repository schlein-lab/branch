#include <gtest/gtest.h>

#include "backend/backend_vtable.hpp"

using branch::backend::Backend;
using branch::backend::BackendVTable;
using branch::backend::BubbleClass;
using branch::backend::ClassificationResult;
using branch::backend::OverlapPair;
using branch::backend::VAFEstimate;

namespace {

// A trivial mock backend used only to exercise the VTable dispatch.
struct MockCtx {
    int destroy_calls = 0;
    int classify_calls = 0;
    int vaf_calls = 0;
    int overlap_calls = 0;
};

void mock_destroy(void* ctx) {
    static_cast<MockCtx*>(ctx)->destroy_calls += 1;
}
void mock_compute_overlaps(void* ctx,
                           const branch::backend::ReadBatch*,
                           std::span<OverlapPair>,
                           std::size_t* out_count) {
    static_cast<MockCtx*>(ctx)->overlap_calls += 1;
    *out_count = 0;
}
void mock_classify_batch(void* ctx,
                         std::span<const branch::backend::BubbleCandidate>,
                         std::span<ClassificationResult> out) {
    static_cast<MockCtx*>(ctx)->classify_calls += 1;
    for (auto& r : out) {
        r.label = BubbleClass::Unknown;
        r.confidence = 0.0F;
    }
}
void mock_estimate_vaf(void* ctx,
                       std::span<const std::uint32_t>,
                       std::span<VAFEstimate> out) {
    static_cast<MockCtx*>(ctx)->vaf_calls += 1;
    for (auto& v : out) {
        v.point = 0.0F;
        v.ci_low = 0.0F;
        v.ci_high = 0.0F;
    }
}
const char* mock_name(void*) { return "mock"; }

constexpr BackendVTable kMockVTable{
    .destroy = &mock_destroy,
    .compute_overlaps = &mock_compute_overlaps,
    .classify_batch = &mock_classify_batch,
    .estimate_vaf_batch = &mock_estimate_vaf,
    .name = &mock_name,
};

}  // namespace

TEST(BackendVTableTest, OverlapPairPacks_to_24_bytes) {
    EXPECT_EQ(sizeof(OverlapPair), 24u);
}

TEST(BackendVTableTest, EmptyBackend_reports_empty) {
    Backend b;
    EXPECT_TRUE(b.empty());
}

TEST(BackendVTableTest, Backend_dispatch_invokes_vtable) {
    MockCtx ctx;
    Backend b(&ctx, &kMockVTable);
    EXPECT_FALSE(b.empty());
    EXPECT_STREQ(b.name(), "mock");

    std::array<ClassificationResult, 2> cls_out{};
    b.classify_batch({}, std::span<ClassificationResult>(cls_out));
    EXPECT_EQ(ctx.classify_calls, 1);
    EXPECT_EQ(cls_out[0].label, BubbleClass::Unknown);

    std::array<VAFEstimate, 1> vaf_out{};
    b.estimate_vaf_batch({}, std::span<VAFEstimate>(vaf_out));
    EXPECT_EQ(ctx.vaf_calls, 1);

    std::size_t n = 0;
    b.compute_overlaps(nullptr, {}, &n);
    EXPECT_EQ(ctx.overlap_calls, 1);
    EXPECT_EQ(n, 0u);
}

TEST(BackendVTableTest, Backend_destructor_calls_destroy_exactly_once) {
    MockCtx ctx;
    {
        Backend b(&ctx, &kMockVTable);
        EXPECT_EQ(ctx.destroy_calls, 0);
    }
    EXPECT_EQ(ctx.destroy_calls, 1);
}

TEST(BackendVTableTest, Move_transfers_ownership_no_double_destroy) {
    MockCtx ctx;
    Backend a(&ctx, &kMockVTable);
    Backend b(std::move(a));
    EXPECT_TRUE(a.empty());
    EXPECT_FALSE(b.empty());
    // dtor of a should not destroy ctx (ownership moved out)
    // dtor of b at end-of-scope destroys ctx exactly once
}
