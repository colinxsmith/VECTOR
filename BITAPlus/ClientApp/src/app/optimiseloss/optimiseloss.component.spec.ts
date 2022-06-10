import { ComponentFixture, TestBed } from '@angular/core/testing';

import { OptimiselossComponent } from './optimiseloss.component';

describe('OptimiselossComponent', () => {
  let component: OptimiselossComponent;
  let fixture: ComponentFixture<OptimiselossComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ OptimiselossComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(OptimiselossComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
